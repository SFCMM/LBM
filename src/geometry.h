#ifndef GRIDGENERATOR_GEOMETRY_H
#define GRIDGENERATOR_GEOMETRY_H

#include <json.h>
#include <memory>
#include <mpi.h>
#include <utility>
#include <vector>
#include <sfcmm_common.h>
#include "common/configuration.h"
#include "functions.h"


template <GInt NDIM>
using Point = VectorD<NDIM>;

using json = nlohmann::json;

// todo: move to boundary condition
enum class BoundaryConditionType { Wall };

class GeometryInterface {
 public:
  GeometryInterface(const MPI_Comm comm) : m_comm(comm){};
  virtual ~GeometryInterface()                                   = default;
  GeometryInterface(const GeometryInterface&)                    = delete;
  GeometryInterface(GeometryInterface&&)                         = delete;
  auto operator=(const GeometryInterface&) -> GeometryInterface& = delete;
  auto operator=(GeometryInterface&&) -> GeometryInterface&      = delete;

  virtual void                      setup(const json& geometry)                                                     = 0;
  virtual inline auto               pointIsInside(const GDouble* x) const -> GBool                                  = 0;
  virtual inline auto               cutWithCell(const GDouble* cellCenter, const GDouble cellLength) const -> GBool = 0;
  [[nodiscard]] virtual inline auto noObjects() const -> GInt                                                       = 0;
  [[nodiscard]] virtual inline auto getBoundingBox() const -> BoundingBoxDynamic                                    = 0;
  [[nodiscard]] virtual auto inline noElements() const -> GInt                                                      = 0;
  [[nodiscard]] virtual auto inline noElements(GInt objId) const -> GInt                                            = 0;

 private:
  MPI_Comm m_comm;
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class GeometryRepresentation {
 public:
  GeometryRepresentation(const json& geom)
    : m_body(config::opt_config_value(geom, "body", GString("unique"))), m_subtract(config::opt_config_value(geom, "subtract", false)){};
  GeometryRepresentation()                                                 = default;
  virtual ~GeometryRepresentation()                                        = default;
  GeometryRepresentation(const GeometryRepresentation&)                    = delete;
  GeometryRepresentation(GeometryRepresentation&&)                         = delete;
  auto operator=(const GeometryRepresentation&) -> GeometryRepresentation& = delete;
  auto operator=(GeometryRepresentation&&) -> GeometryRepresentation&      = delete;

  [[nodiscard]] virtual inline auto pointIsInside(const Point<NDIM>& x) const -> GBool                        = 0;
  [[nodiscard]] virtual inline auto cutWithCell(const Point<NDIM>& center, GDouble cellLength) const -> GBool = 0;
  [[nodiscard]] virtual inline auto getBoundingBox() const -> BoundingBoxDynamic                              = 0;
  [[nodiscard]] virtual inline auto str() const -> GString                                                    = 0;
  [[nodiscard]] virtual inline auto noElements() const -> GInt                                                = 0;
  [[nodiscard]] virtual inline auto min(const GInt dir) const -> GDouble                                      = 0;
  [[nodiscard]] virtual inline auto max(const GInt dir) const -> GDouble                                      = 0;

  [[nodiscard]] inline auto type() const -> GeomType { return m_type; }
  // necessary if objectref cannot be cast to const
  [[nodiscard]] inline auto ctype() const -> GeomType { return m_type; }
  [[nodiscard]] inline auto name() const -> GString { return m_name; }
  // necessary if objectref cannot be cast to const
  [[nodiscard]] inline auto cname() const -> GString { return m_name; }
  [[nodiscard]] inline auto subtract() const -> GBool { return m_subtract; }
  [[nodiscard]] inline auto body() const -> GString { return m_body; }
  inline auto               body() -> GString& { return m_body; }
  [[nodiscard]] inline auto elementOffset() const -> GInt { return m_elementOffset; }
  inline auto               elementOffset() -> GInt& { return m_elementOffset; }

 protected:
  inline auto type() -> GeomType& { return m_type; }
  inline auto name() -> GString& { return m_name; }

 private:
  GeomType              m_type = GeomType::unknown;
  GString               m_name = "undefined";
  GString               m_body;
  GBool                 m_subtract          = false;
  BoundaryConditionType m_boundaryCondition = BoundaryConditionType::Wall;
  GInt                  m_elementOffset     = -1;
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class GeometrySTL : public GeometryRepresentation<DEBUG_LEVEL, NDIM> {
 public:
  GeometrySTL(GString fileName, const GString& _name) : m_fileName(std::move(fileName)) {
    name() = _name;
    type() = GeomType::stl;
    loadFile();
  };

  GeometrySTL(const json& stl, const GString& _name) : GeometryRepresentation<DEBUG_LEVEL, NDIM>(stl), m_fileName(stl["filename"]) {
    name() = _name;
    type() = GeomType::stl;
    loadFile();
  };

  [[nodiscard]] auto inline pointIsInside(const Point<NDIM>& x) const -> GBool override {
    // 1. check if it is in object bounding box
    if(!pointInsideObjBB(x)) {
      return false;
    }


    // 2. cast ray from the point to the outside and determine cuts (if an uneven number of cuts is found the point is inside)
    std::vector<GInt>             nodeList;
    std::array<GDouble, 2 * NDIM> targetRegion;
    for(GInt dir = 0; dir < NDIM; ++dir) {
      targetRegion[2 * dir]     = x[dir];
      targetRegion[2 * dir + 1] = x[dir];
    }


    static constexpr GDouble tolerance = 1E-10;
    std::vector<Point<NDIM>> intersectionPoints;

    for(GInt dir = 0; dir < NDIM; ++dir) {
      // cast ray to the outside (too make sure 2*the extend)
      targetRegion[2 * dir + 1] += 2 * m_extend[dir];
      m_kd.retrieveNodes(targetRegion, nodeList);
      // reset
      targetRegion[2 * dir + 1] = x[dir];

      // it is enough to cast one ray in any direction otherwise the geometry has holes!
      Point<NDIM> ray;
      ray.fill(0);
      ray[dir] = 2 * m_extend[dir];

      for(const GInt triId : nodeList) {
        const auto& tri   = m_triangles[triId];
        const auto  edge3 = x - tri.m_vertices[0];


        const auto a = -tri.m_normal.dot(edge3);
        const auto b = tri.m_normal.dot(ray);

        if(fabs(b) < tolerance) {
          if(fabs(a) < tolerance) {
            // ray lies in triangle plane
            logger << "Found ray in triangle plane" << std::endl;
            return true;
          }
          // ray disjoint from plane
          continue;
        }
        const auto r = a / b;

        // check necessary conditions
        if(r < -tolerance || r > 1 + tolerance) {
          continue;
        }

        const auto ip = x + r * ray;
        const auto w  = ip - tri.m_vertices[0];
        const auto u  = tri.m_vertices[1] - tri.m_vertices[0];
        const auto uu = u.dot(u);
        const auto wu = w.dot(u);
        const auto v  = tri.m_vertices[2] - tri.m_vertices[0];
        const auto vv = v.dot(v);
        const auto uv = u.dot(v);
        const auto wv = w.dot(v);
        const auto D  = uv * uv - uu * vv;
        const auto s  = (uv * wv - vv * wu) / D;
        if(s < -tolerance || s > 1 + tolerance) {
          // Intersection point is outside triangle
          continue;
        }

        const auto t = (uv * wu - uu * wv) / D;
        if(t < -tolerance || (s + t) > 1 + tolerance) {
          // Intersection point is outside triangle
          continue;
        }

        if(intersectionPoints.empty()) {
          intersectionPoints.emplace_back(ip);
        } else {
          GBool singleIP = true;

          // check if the intersection point already exists (i.e. shared edges)
          for(const auto& ip2 : intersectionPoints) {
            const auto diff = ip - ip2;
            if(diff.norm() < tolerance) {
              singleIP = false;
              break;
            }
          }
          if(singleIP) {
            intersectionPoints.emplace_back(ip);
          }
        }
      }

      if(isEven(intersectionPoints.size())) {
        return false;
      }

      intersectionPoints.clear();
      nodeList.clear();
    }

    return true;
  }

  [[nodiscard]] inline auto cutWithCell(const Point<NDIM>& cellCenter, const GDouble cellLength) const -> GBool override {
    if(cellCutWithObjBB(cellCenter, cellLength)) {
      std::vector<GInt>             nodeList;
      std::array<GDouble, 2 * NDIM> targetRegion;
      for(GInt dir = 0; dir < NDIM; ++dir) {
        // search for cuts within the bb of the current cell
        targetRegion[2 * dir]     = cellCenter[dir] - HALF * cellLength;
        targetRegion[2 * dir + 1] = cellCenter[dir] + HALF * cellLength;
      }

      // obtain kd tree nodes which have possible cuts
      m_kd.retrieveNodes(targetRegion, nodeList);
      if(DEBUG_LEVEL > Debug_Level::debug) {
        if(!nodeList.empty()) {
          logger << "possible nodes " << strStreamify(nodeList).str() << std::endl;
        } else {
          logger << "no nodes!" << std::endl;
        }
      }
      removeNonOverlappingNodes(targetRegion, nodeList);

      if(DEBUG_LEVEL > Debug_Level::debug) {
        if(!nodeList.empty()) {
          logger << "possible nodes after non-overlapping removal " << strStreamify(nodeList).str() << std::endl;
        } else {
          logger << "no nodes left after non-overlapping removal!" << std::endl;
        }
      }

      if(!nodeList.empty()) {
        static constexpr GInt combo[3][6]    = {{0, 2, 0, 2, 1, 2}, {0, 2, 0, 2, 0, 1}, {0, 1, 0, 1, 1, 2}};
        const GDouble         cellHalfLength = HALF * cellLength;
        for(const GInt triId : nodeList) {
          GBool                            cut  = true;
          auto&                            tri  = m_triangles[triId];
          const std::array<Point<NDIM>, 3> vert = {tri.m_vertices[0] - cellCenter, tri.m_vertices[1] - cellCenter,
                                                   tri.m_vertices[2] - cellCenter};
          //          const Point<NDIM>                res  = vert[1] - vert[0];
          const std::array<Point<NDIM>, 3> edge = {vert[1] - vert[0], vert[2] - vert[1], vert[0] - vert[2]};


          for(GInt i = 0; i < 3; ++i) {
            GInt mul = 1;
            GInt j   = 2;
            GInt k   = 1;
            for(GInt l = 0; l < NDIM; ++l) {
              GDouble p0 = mul * edge[i][j] * vert[combo[i][2 * l]][k] - mul * edge[i][k] * vert[combo[i][2 * l]][j];
              GDouble p1 = mul * edge[i][j] * vert[combo[i][2 * l + 1]][k] - mul * edge[i][k] * vert[combo[i][2 * l + 1]][j];

              const GDouble min = (p0 < p1) ? p0 : p1;
              const GDouble max = (p0 < p1) ? p1 : p0;
              const GDouble rad = cellHalfLength * (abs(edge[i][j]) + abs(edge[i][k]));

              if(min > rad || max < -rad) {
                cut = false;
                break;
              }

              mul *= -1;
              j -= l;
              k = 0;
            }
            if(!cut) {
              break;
            }
          }
          if(!cut) {
            continue;
          }

          for(GInt i = 0; i < NDIM; i++) {
            GDouble min = vert[0][i];
            GDouble max = vert[0][i];

            for(GInt j = 1; j < NDIM; j++) {
              if(vert[j][i] < min) {
                min = vert[j][i];
              }
              if(vert[j][i] > max) {
                max = vert[j][i];
              }
            }

            if(min > cellHalfLength || max < -cellHalfLength) {
              cut = false;
              break;
            }
          }
          if(!cut) {
            continue;
          }
          Point<NDIM> normal;
          // todo: fix for 1D/2D
          normal[0] = edge[0][1] * edge[1][2] - edge[0][2] * edge[1][1];
          normal[1] = edge[0][2] * edge[1][0] - edge[0][0] * edge[1][2];
          normal[2] = edge[0][0] * edge[1][1] - edge[0][1] * edge[1][0];

          GDouble d = normal.dot(vert[0]);
          d *= -1;

          Point<NDIM> vmin;
          Point<NDIM> vmax;
          for(GInt i = 0; i < NDIM; i++) {
            if(normal[i] > 0.0) {
              vmin[i] = -cellHalfLength;
              vmax[i] = cellHalfLength;
            } else {
              vmin[i] = cellHalfLength;
              vmax[i] = -cellHalfLength;
            }
          }

          if((normal[0] * vmin[0] + normal[1] * vmin[1] + normal[2] * vmin[2]) + d > 0.0) {
            continue;
          }
          if((normal[0] * vmax[0] + normal[1] * vmax[1] + normal[2] * vmax[2]) + d >= 0.0) {
            return true;
          }
        }
      }
    }

    return false;
  }

  [[nodiscard]] inline auto getBoundingBox() const -> BoundingBoxDynamic override { return BoundingBoxDynamic(m_bbox); }

  [[nodiscard]] inline auto pointInsideObjBB(const Point<NDIM>& x) const -> GBool {
    for(GInt dir = 0; dir < NDIM; ++dir) {
      const GDouble p = x[dir];
      if(p < m_bbox.min(dir) || p > m_bbox.max(dir)) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] inline auto cellCutWithObjBB(const Point<NDIM>& cellCenter, const GDouble cellLength) const -> GBool {
    for(GInt dir = 0; dir < NDIM; ++dir) {
      const GDouble diffMin = abs(cellCenter[dir] - m_bbox.min(dir));
      const GDouble diffMax = abs(cellCenter[dir] - m_bbox.max(dir));
      //      cerr0 << "diffMin " << diffMin << " diffMax " << diffMax << std::endl;
      if(diffMin <= cellLength || diffMax <= cellLength) {
        // is closer than cellLength a cut might exist
        return true;
      }
    }
    // point is within overall BB so might cut the geometry
    return pointInsideObjBB(cellCenter);
  }

  [[nodiscard]] inline auto noElements() const -> GInt override { return m_noTriangles; }

  [[nodiscard]] inline auto str() const -> GString override {
    std::stringstream ss;
    ss << SP1 << "STL"
       << "\n";
    ss << SP7 << "Name: " << name() << "\n";
    ss << SP7 << "Body: " << body() << "\n";
    ss << SP7 << "Filename: " << m_fileName << "\n";
    ss << SP7 << "Binary: " << std::boolalpha << m_binary << "\n";
    ss << SP7 << "No triangles: " << m_noTriangles << "\n";
    ss << SP7 << "Bounding Box: " << m_bbox.str() << "\n";
    ss << SP7 << "Extend: " << strStreamify<NDIM>(m_extend).str() << "\n";
    return ss.str();
  }

  [[nodiscard]] inline auto min(const GInt dir) const -> GDouble override { return m_bbox.min(dir); }
  [[nodiscard]] inline auto max(const GInt dir) const -> GDouble override { return m_bbox.max(dir); }

  void printElements() const {
    GInt elementId = 0;
    for(const auto& tri : m_triangles) {
      std::cout << "----- Element " << elementId++ << "-----" << std::endl;
      triangle_::print(tri);
    }
  }

  inline auto triangles() const -> const std::vector<triangle<NDIM>>& { return m_triangles; }

 private:
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::name;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::body;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::type;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::elementOffset;

  void loadFile() {
    checkFileExistence();
    determineBinary();
    countTriangles();
    m_triangles.resize(m_noTriangles);
    if(m_binary) {
      readBinarySTL();
    } else {
      readASCIISTL();
    }
    determineBoundaryBox();
    m_kd.buildTree(m_triangles, m_bbox);
  }

  void checkFileExistence() {
    if(!isFile(m_fileName)) {
      TERMM(-1, "The STL file: " + m_fileName + " cannot be found!");
    }
  }

  void determineBinary() {
    // A STL file is ASCII if the first 5 letters in a file are "solid"!
    std::ifstream                  ifl(m_fileName);
    static constexpr GInt          buffer_size = 1024;
    std::array<GChar, buffer_size> buffer{};
    ifl.getline(buffer.data(), buffer_size);
    if(static_cast<GString>(buffer.data()).find("solid") != 0) {
      m_binary = true;
    }
    ifl.close();
  }

  void countTriangles() {
    if(m_binary) {
      static constexpr GInt stl_header_size   = 80 + 4;
      static constexpr GInt stl_triangle_size = 50;
      m_noTriangles                           = (fileSize(m_fileName) - stl_header_size) / stl_triangle_size;
    } else {
      // open file in ascii format
      std::ifstream ifl(m_fileName);

      static constexpr GInt          buffer_size = 1024;
      std::array<GChar, buffer_size> buffer{};

      GString text;
      // If number of elements unknown, count them...
      while(text.find("endsolid") == std::string::npos) {
        ifl.getline(buffer.data(), buffer_size);
        text = static_cast<GString>(buffer.data());
        if(text.find("facet") != std::string::npos && text.find("endface") == std::string::npos) {
          ++m_noTriangles;
        }
      }
      ifl.close();
    }
  }

  void readASCIISTL() {
    // open file in ascii format
    std::ifstream ifl(m_fileName);

    static constexpr GInt          buffer_size = 1024;
    std::array<GChar, buffer_size> buffer{};
    GString                        text;

    // iterate overall triangles
    GInt triangleId = 0;
    while(text.find("endsolid") == std::string::npos) {
      ifl.getline(buffer.data(), buffer_size);
      text = static_cast<GString>(buffer.data());

      // get normal of the triangle
      if(text.find("normal") != std::string::npos) {
        trim(text);

        // split the string into multiple delimited by " "
        // todo: note this is not actually correct STL value can be split also by tabs!
        std::vector<std::string> tokens;
        tokenize(text, tokens, " ", true);

        // Read normal vector components
        for(GInt i = 0; i < NDIM; i++) {
          tokens[i + 2] = trim(tokens[i + 2]);
          // todo: role into function
          if(tokens[i + 2].find_first_not_of("0123456789.eE-+") != std::string::npos) {
            TERMM(-1, "ERROR: normal component " + tokens[i + 2] + " is not a number ");
          }

          m_triangles[triangleId].m_normal[i] = stod(tokens[i + 2]);
        }
        // Jump expression : outer loop
        ifl.getline(buffer.data(), buffer_size);

        // Read vertices
        for(GInt j = 0; j < NDIM; j++) {
          ifl.getline(buffer.data(), buffer_size);

          text = static_cast<GString>(buffer.data());
          trim(text);
          std::vector<std::string> vertex;
          tokenize(text, vertex, " ", true);

          for(GInt i = 0; i < NDIM; i++) {
            vertex[i + 1] = trim(vertex[i + 1]);
            // todo: role into function
            if(vertex[i + 1].find_first_not_of("0123456789.eE-+") != std::string::npos) {
              TERMM(-1, "ERROR: vertex component " + vertex[i + 1] + " is not a number ");
            }


            m_triangles[triangleId].m_vertices[j][i] = stod(vertex[i + 1]);
          }
        }

        setMinMax(m_triangles[triangleId]);
        ++triangleId;
      }
    }
    ifl.close();
  }

  void readBinarySTL() {
    using namespace std;

    // data in an binary STL file normal, vertex1, vertex2, vertex3
    using facet_t = struct {
      std::array<float, NDIM> n, v1, v2, v3;
    };
    facet_t facet;

    // open file in binary format
    FILE* fp = fopen(m_fileName.c_str(), "rb");

    static constexpr GInt        stl_header_size = 85;
    array<char, stl_header_size> header{};
    unsigned short               ibuff2 = 0;

    if(fread(header.data(), 1, stl_header_size - 1, fp) == 0U) {
      TERMM(-1, "ERROR: Memory error!");
    }

    GInt                  triangleId      = 0;
    static constexpr GInt float_byte_size = 4;
    static constexpr GInt no_values       = 4; // (3 vertices + 1 normal)
    for(GInt i = 0; fread(&facet, no_values * float_byte_size * NDIM, 1, fp) > 0; i++) {
      if(fread(&ibuff2, 2, 1, fp) == 0U) {
        TERMM(-1, "ERROR: Memory error!");
      }

      for(GInt j = 0; j < NDIM; j++) {
        m_triangles[triangleId].m_normal[j]      = facet.n[j];
        m_triangles[triangleId].m_vertices[0][j] = facet.v1[j];
        m_triangles[triangleId].m_vertices[1][j] = facet.v2[j];
        m_triangles[triangleId].m_vertices[2][j] = facet.v3[j];
      }

      setMinMax(m_triangles[triangleId]);
      ++triangleId;
    }

    fclose(fp);
  }

  void setMinMax(triangle<NDIM>& tri) {
    for(GInt j = 0; j < NDIM; j++) {
      tri.m_max[j] = tri.m_vertices[0][j];
      tri.m_min[j] = tri.m_vertices[0][j];
    }

    for(GInt vertexId = 0; vertexId < 3; vertexId++) {
      for(GInt dir = 0; dir < NDIM; dir++) {
        const GDouble value = tri.m_vertices[vertexId][dir];
        // Find maximum
        tri.m_max[dir] = (tri.m_max[dir] < value) ? value : tri.m_max[dir];
        // Find minimum
        tri.m_min[dir] = (tri.m_min[dir] > value) ? value : tri.m_min[dir];
      }
    }
    for(GInt dir = 0; dir < NDIM; dir++) {
      tri.m_max[dir] += GDoubleEps;
      tri.m_min[dir] -= GDoubleEps;
    }
  }

  void determineBoundaryBox() {
    // initialize
    for(GInt dir = 0; dir < NDIM; ++dir) {
      m_bbox.min(dir) = m_triangles[0].m_min[dir];
      m_bbox.max(dir) = m_triangles[0].m_max[dir];
    }

    for(const auto& tri : m_triangles) {
      for(GInt dir = 0; dir < NDIM; dir++) {
        // Find minimum
        m_bbox.min(dir) = (m_bbox.min(dir) > tri.m_min[dir]) ? tri.m_min[dir] : m_bbox.min(dir);
        // Find maximum
        m_bbox.max(dir) = (m_bbox.max(dir) < tri.m_max[dir]) ? tri.m_max[dir] : m_bbox.max(dir);
      }
    }

    for(GInt dir = 0; dir < NDIM; dir++) {
      m_extend[dir] = m_bbox.max(dir) - m_bbox.min(dir);
    }
  }

  void removeNonOverlappingNodes(const std::array<GDouble, 2 * NDIM>& targetRegion, std::vector<GInt>& nodeList) const {
    for(auto it = nodeList.begin(); it != nodeList.end();) {
      const auto& tri = m_triangles[*it];

      GBool erase = false;
      for(GInt i = 0; i < NDIM; i++) {
        //        cerr0 << "tri " << tri.m_min[i] << " > " << targetRegion [2 * i + 1] << " " << tri.m_max[i] << " < " << targetRegion[2*i]
        //        <<
        //            std::endl;
        if(tri.m_min[i] > targetRegion[2 * i + 1] || tri.m_max[i] < targetRegion[2 * i]) {
          //          cerr0 << "erase " << std::endl;
          it    = nodeList.erase(it);
          erase = true;
          break;
        }
      }
      if(!erase) {
        ++it;
      }
    }
  }

  GString                     m_fileName;
  GBool                       m_binary      = false; // file is ASCII or binary
  GInt                        m_noTriangles = 0;
  std::vector<triangle<NDIM>> m_triangles;
  BoundingBoxCT<NDIM>         m_bbox;
  std::array<GDouble, NDIM>   m_extend{};
  KDTree<DEBUG_LEVEL, NDIM>   m_kd;
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class GeometryAnalytical : public GeometryRepresentation<DEBUG_LEVEL, NDIM> {
 public:
  GeometryAnalytical() = default;
  GeometryAnalytical(const json& geom) : GeometryRepresentation<DEBUG_LEVEL, NDIM>(geom){};

  [[nodiscard]] inline auto noElements() const -> GInt override { return 1; }
  [[nodiscard]] inline auto min(const GInt dir) const -> GDouble override { return this->getBoundingBox().min(dir); }
  [[nodiscard]] inline auto max(const GInt dir) const -> GDouble override { return this->getBoundingBox().max(dir); }

 private:
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class GeomSphere : public GeometryAnalytical<DEBUG_LEVEL, NDIM> {
 public:
  GeomSphere(const Point<NDIM>& center, GDouble radius, const GString& _name) : m_center(center), m_radius(radius) {
    name() = _name;
    type() = GeomType::sphere;
  };
  GeomSphere(const json& sphere, const GString& _name)
    : GeometryAnalytical<DEBUG_LEVEL, NDIM>(sphere),
      m_center(static_cast<std::vector<GDouble>>(sphere["center"]).data()),
      m_radius(sphere["radius"]) {
    name() = _name;
    type() = GeomType::sphere;
  };

  [[nodiscard]] auto inline pointIsInside(const Point<NDIM>& x) const -> GBool override {
    return (x - m_center).norm() < m_radius + GDoubleEps;
  }

  [[nodiscard]] inline auto cutWithCell(const Point<NDIM>& cellCenter, GDouble cellLength) const -> GBool override {
    const GDouble distance        = (cellCenter - m_center).norm();
    const GDouble halfCellLength  = 0.5 * cellLength;
    const GDouble surrCircleCellR = gcem::sqrt(2) * halfCellLength;

    // cell has a possible cut with the sphere
    if(distance <= m_radius + surrCircleCellR && distance >= m_radius - surrCircleCellR) {
      GInt          pointsInside = 0;
      VectorD<NDIM> l;
      l.fill(halfCellLength);
      Point<NDIM> cellOriginVert = cellCenter - l;

      if constexpr(NDIM == 2) {
        VectorD<2> ex = {1, 0};
        VectorD<2> ey = {0, 1};

        for(GInt dirX = 0; dirX < NDIM; ++dirX) {
          for(GInt dirY = 0; dirY < NDIM; ++dirY) {
            const Point<NDIM> vert  = cellOriginVert + cellLength * dirX * ex + cellLength * dirY * ey;
            const GDouble     vertD = (vert - m_center).norm();
            if(vertD <= m_radius) {
              ++pointsInside;
            }
          }
        }
      }
      if constexpr(NDIM == 3) {
        TERMM(-1, "not implemented");
      }
      if(pointsInside < gcem::pow(2, NDIM)) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] inline auto getBoundingBox() const -> BoundingBoxDynamic override {
    BoundingBoxDynamic bbox;
    bbox.init(NDIM);

    for(GInt dir = 0; dir < NDIM; ++dir) {
      bbox.min(dir) = m_center[dir] - m_radius;
      bbox.max(dir) = m_center[dir] + m_radius;
    }
    return bbox;
  }

  [[nodiscard]] inline auto str() const -> GString override {
    std::stringstream ss;
    ss << SP1 << "Sphere"
       << "\n";
    ss << SP7 << "Name: " << name() << "\n";
    ss << SP7 << "Body: " << body() << "\n";
    ss << SP7 << "Center: " << strStreamify<NDIM>(m_center).str() << "\n";
    ss << SP7 << "Radius: " << m_radius << "\n";
    return ss.str();
  }

 private:
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::name;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::body;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::type;

  Point<NDIM> m_center{NAN};
  GDouble     m_radius = 0;
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class GeomBox : public GeometryAnalytical<DEBUG_LEVEL, NDIM> {
 public:
  GeomBox(const Point<NDIM>& A, const Point<NDIM>& B, const GString& _name) : m_A(A), m_B(B) {
    name() = _name;
    type() = GeomType::box;
    checkValid();
  }
  GeomBox(const json& box, const GString& _name) : GeometryAnalytical<DEBUG_LEVEL, NDIM>(box) {
    const auto A = std::vector<GDouble>(box["A"]);
    const auto B = std::vector<GDouble>(box["B"]);

    if(A.size() != NDIM || B.size() != NDIM) {
      TERMM(-1, "Invalid dimensionality given for box corner points");
    }

    m_A = Point<NDIM>(A.data());
    m_B = Point<NDIM>(B.data());

    name() = _name;
    type() = GeomType::box;
    checkValid();
  }

  [[nodiscard]] auto inline pointIsInside(const Point<NDIM>& x) const -> GBool override {
    for(GInt dir = 0; dir < NDIM; ++dir) {
      if(m_A[dir] > x[dir] || m_B[dir] < x[dir]) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] inline auto cutWithCell(const Point<NDIM>& cellCenter, GDouble cellLength) const -> GBool override {
    const GDouble halfCellLength = 0.5 * cellLength;

    if constexpr(NDIM == 1) {
      if(std::abs(m_A[0] - cellCenter[0]) <= halfCellLength) {
        return true;
      }
      if(std::abs(m_B[0] - cellCenter[0]) <= halfCellLength) {
        return true;
      }
    }


    if constexpr(NDIM == 2) {
      if(m_A[0] - halfCellLength <= cellCenter[0] and m_B[0] + halfCellLength >= cellCenter[0]) {
        if(std::abs(cellCenter[1] - m_A[1]) <= halfCellLength) {
          return true;
        }
        if(std::abs(cellCenter[1] - m_B[1]) <= halfCellLength) {
          return true;
        }
      }
      if(m_A[1] - halfCellLength <= cellCenter[1] and m_B[1] + halfCellLength >= cellCenter[1]) {
        if(std::abs(cellCenter[0] - m_A[0]) <= halfCellLength) {
          return true;
        }
        if(std::abs(cellCenter[0] - m_B[0]) <= halfCellLength) {
          return true;
        }
      }
    }

    if constexpr(NDIM == 3) {
      TERMM(-1, "impl");
    }


    //    for(GInt dir = 0; dir < NDIM; ++dir) {
    //      if((m_A[dir] > cellCenter[dir] || m_B[dir] < cellCenter[dir]) // cellCenter is within the box
    //         && abs(m_A[dir] - cellCenter[dir]) > cellLength && abs(m_B[dir] - cellCenter[dir]) > cellLength /*cellcenter cuts the box*/)
    //         {
    //        return false;
    //      }
    //    }
    return false;
  }

  [[nodiscard]] inline auto getBoundingBox() const -> BoundingBoxDynamic override {
    BoundingBoxDynamic bbox;
    bbox.init(NDIM);

    for(GInt dir = 0; dir < NDIM; ++dir) {
      bbox.min(dir) = m_A[dir];
      bbox.max(dir) = m_B[dir];
    }
    return bbox;
  }

  [[nodiscard]] inline auto str() const -> GString override {
    std::stringstream ss;
    ss << SP1 << "Box"
       << "\n";
    ss << SP7 << "Name: " << name() << "\n";
    ss << SP7 << "Body: " << body() << "\n";
    ss << SP7 << "Point A: " << strStreamify<NDIM>(m_A).str() << "\n";
    ss << SP7 << "Point B: " << strStreamify<NDIM>(m_B).str() << "\n";
    ss << SP7 << "Subtract: " << std::boolalpha << subtract() << "\n";
    return ss.str();
  }

 private:
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::name;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::body;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::type;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::subtract;

  void checkValid() const {
    for(GInt dir = 0; dir < NDIM; dir++) {
      if(m_A[dir] > m_B[dir]) {
        TERMM(-1, "ERROR: The specification of the box is invalid " + std::to_string(m_A[dir]) + " > " + std::to_string(m_B[dir]));
      }
    }
  }
  Point<NDIM> m_A;
  Point<NDIM> m_B;
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class GeomCube : public GeometryAnalytical<DEBUG_LEVEL, NDIM> {
 public:
  GeomCube(const Point<NDIM>& center, const GDouble length, const GString& _name) : m_center(center), m_length(length) {
    name() = _name;
    type() = GeomType::cube;
  };
  GeomCube(const json& cube, const GString& _name)
    : GeometryAnalytical<DEBUG_LEVEL, NDIM>(cube),
      m_center(static_cast<std::vector<GDouble>>(cube["center"]).data()),
      m_length(cube["length"]) {
    name() = _name;
    type() = GeomType::cube;
  };

  [[nodiscard]] auto inline pointIsInside(const Point<NDIM>& x) const -> GBool override {
    for(GInt dir = 0; dir < NDIM; ++dir) {
      if(abs(x[dir] - m_center[dir]) > m_length) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] inline auto cutWithCell(const Point<NDIM>& cellCenter, GDouble cellLength) const -> GBool override {
    for(GInt dir = 0; dir < NDIM; ++dir) {
      if(abs(cellCenter[dir] - m_center[dir]) > m_length + cellLength) {
        return false;
      }
    }
    return true;
  }

  [[nodiscard]] inline auto getBoundingBox() const -> BoundingBoxDynamic override {
    BoundingBoxDynamic bbox;
    bbox.init(NDIM);

    GDouble cicumference_radius = gcem::sqrt(NDIM) * m_length;
    for(GInt dir = 0; dir < NDIM; ++dir) {
      bbox.min(dir) = m_center[dir] - cicumference_radius;
      bbox.max(dir) = m_center[dir] + cicumference_radius;
    }
    return bbox;
  }

  [[nodiscard]] inline auto str() const -> GString override {
    std::stringstream ss;
    ss << SP1 << "Cube"
       << "\n";
    ss << SP7 << "Name: " << name() << "\n";
    ss << SP7 << "Body: " << body() << "\n";
    ss << SP7 << "Center: " << strStreamify<NDIM>(m_center).str() << "\n";
    ss << SP7 << "Length: " << m_length << "\n";
    return ss.str();
  }


 private:
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::name;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::body;
  using GeometryRepresentation<DEBUG_LEVEL, NDIM>::type;

  Point<NDIM> m_center{NAN};
  GDouble     m_length = 0;
};

struct MinMaxType {
  std::function<GDouble(GInt)> min;
  std::function<GDouble(GInt)> max;
};

template <Debug_Level DEBUG_LEVEL, GInt NDIM>
class GeometryManager : public GeometryInterface {
 public:
  GeometryManager(const MPI_Comm comm) : GeometryInterface(comm){};

  void setup(const json& geometry) override {
    // generate geometric representation objects for each defined geometry
    for(const auto& object : geometry.items()) {
      const GString& name = object.key();

      if(name == "0") {
        TERMM(-1, "Malformed json: No valid key!");
      }

      if(geometry[name].count("type") == 0) {
        TERMM(-1, "Malformed json: No \"type\" given!");
      }

      switch(resolveGeomType(geometry[name]["type"])) {
        case GeomType::sphere: {
          m_geomObj.emplace_back(std::make_unique<GeomSphere<DEBUG_LEVEL, NDIM>>(geometry[name], name));
          break;
        }
        case GeomType::box: {
          m_geomObj.emplace_back(std::make_unique<GeomBox<DEBUG_LEVEL, NDIM>>(geometry[name], name));
          break;
        }
        case GeomType::cube: {
          m_geomObj.emplace_back(std::make_unique<GeomCube<DEBUG_LEVEL, NDIM>>(geometry[name], name));
          break;
        }
        case GeomType::stl: {
          m_geomObj.template emplace_back(std::make_unique<GeometrySTL<DEBUG_LEVEL, NDIM>>(geometry[name], name));
          break;
        }
        case GeomType::unknown:
          [[fallthrough]];
        default: {
          logger << SP2 << "Unknown geometry type " << object.key() << std::endl;
          break;
        }
      }
    }

    // map the geometry objects to bodies
    for(const auto& geom : m_geomObj) {
      const GString name = geom->cname();
      if(geom->body() == "unique") {
        geom->body() = name;
        m_bodyMap.emplace(name, std::vector<GeometryRepresentation<DEBUG_LEVEL, NDIM>*>({geom.get()}));
      } else {
        const GString body = geom->body();
        if(m_bodyMap.count(body) == 1) {
          m_bodyMap[body].emplace_back(geom.get());
        } else {
          m_bodyMap.emplace(body, std::vector<GeometryRepresentation<DEBUG_LEVEL, NDIM>*>({geom.get()}));
        }
      }
      logger << geom->str() << std::endl;
    }

    // build elementMinMax
    GInt offsetCounter = 0;
    for(auto& geom : m_geomObj) {
      if(geom->ctype() == GeomType::stl) {
        for(const auto& tri : static_cast<GeometrySTL<DEBUG_LEVEL, NDIM>*>(geom.get())->triangles()) {
          std::function<GDouble(GInt)> triMin = [&](GInt dir) { return triangle_::min(tri, dir); };
          std::function<GDouble(GInt)> triMax = [&](GInt dir) { return triangle_::max(tri, dir); };

          m_elementMinMax.emplace_back(MinMaxType{triMin, triMax});
        }
      } else {
        std::function<GDouble(GInt)> geomMin = [&](GInt dir) { return geom->min(dir); };
        std::function<GDouble(GInt)> geomMax = [&](GInt dir) { return geom->max(dir); };

        m_elementMinMax.emplace_back(MinMaxType{geomMin, geomMax});
      }
      geom->elementOffset() = offsetCounter;
      offsetCounter += geom->noElements();
    }

    m_kd.buildTree(*this);
  }

  [[nodiscard]] auto inline pointIsInside(const GDouble* x) const -> GBool override {
    Point<NDIM> point(x);
    return pointIsInside(point);
  }

  [[nodiscard]] auto inline pointIsInside(const Point<NDIM>& point) const -> GBool {
    // \todo: check overall bounding box first

    // iterate through bodies
    for(const auto& [key, geomObjList] : m_bodyMap) {
      // todo: move to global data to not always check
      GBool hasSubtraction = false;
      for(const auto& obj : geomObjList) {
        if(obj->subtract()) {
          hasSubtraction = true;
        }
      }
      // check geometries within body
      for(const auto& obj : geomObjList) {
        if(obj->pointIsInside(point)) {
          // no subtraction geometries in this body so no further checks
          if(!hasSubtraction) {
            return true;
          }

          // geometry should be subtracted
          if(obj->subtract()) {
            return false;
          }

          for(const auto& sub_obj : geomObjList) {
            // point is within subtracted volume so not inside!
            if(sub_obj->subtract() && sub_obj->pointIsInside(point)) {
              return false;
            }
          }

          return true;
        }
      }
    }
    return false;
  }

  [[nodiscard]] auto inline cutWithCell(const GDouble* cellCenter, const GDouble cellLength) const -> GBool override {
    Point<NDIM> center(cellCenter);
    return cutWithCell(center, cellLength);
  }

  [[nodiscard]] auto inline cutWithCell(const Point<NDIM>& cellCenter, const GDouble cellLength) const -> GBool {
    // \todo: check overall bounding box first
    for(const auto& obj : m_geomObj) {
      if(obj->cutWithCell(cellCenter, cellLength)) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] auto inline cutWithCell(const GString& geomName, const Point<NDIM>& cellCenter, const GDouble cellLength) const -> GBool {
    // \todo: check overall bounding box first
    for(const auto& obj : m_geomObj) {
      if(obj->cname() == geomName) {
        if(obj->cutWithCell(cellCenter, cellLength)) {
          return true;
        }
      }
    }
    return false;
  }

  [[nodiscard]] auto inline noObjects() const -> GInt override { return m_geomObj.size(); }

  [[nodiscard]] auto inline noElements() const -> GInt override {
    return std::transform_reduce(
        m_geomObj.begin(), m_geomObj.end(), 0, [](GInt a, GInt b) { return a + b; }, [](const auto& obj) { return obj->noElements(); });
  }

  [[nodiscard]] auto inline noElements(GInt objId) const -> GInt override { return m_geomObj[objId]->noElements(); }

  [[nodiscard]] inline auto getBoundingBox() const -> BoundingBoxDynamic override {
    BoundingBoxDynamic bbox;
    bbox.init(NDIM);

    for(GInt objId = 0; objId < static_cast<GInt>(m_geomObj.size()); ++objId) {
      const auto& obj       = m_geomObj.at(objId);
      const auto  temp_bbox = obj->getBoundingBox();
      for(GInt dir = 0; dir < NDIM; ++dir) {
        // set the bounding box to the values of the first object for initialization
        if(bbox.min(dir) > temp_bbox.min(dir) || objId == 0) {
          bbox.min(dir) = temp_bbox.min(dir);
        }
        if(bbox.max(dir) < temp_bbox.max(dir) || objId == 0) {
          bbox.max(dir) = temp_bbox.max(dir);
        }
      }
    }
    return bbox;
  }

  [[nodiscard]] inline auto elementBoundingBox(const GInt elementId, const GInt dir) const -> GDouble {
    ASSERT(dir <= 2 * NDIM, "Invalid dir");
    if(dir >= NDIM) {
      return m_elementMinMax.at(elementId).max(dir - NDIM);
    }
    return m_elementMinMax.at(elementId).min(dir);
  }

  [[nodiscard]] inline auto pointInElementBB(const GInt elementId, const Point<NDIM>& x) const -> GBool {
    for(GInt dir = 0; dir < NDIM; ++dir) {
      if(x[dir] > m_elementMinMax[elementId].max(dir) || x[dir] < m_elementMinMax[elementId].min(dir)) {
        return false;
      }
    }
    return true;
  }


 private:
  std::vector<std::unique_ptr<GeometryRepresentation<DEBUG_LEVEL, NDIM>>>              m_geomObj;
  std::unordered_map<GString, std::vector<GeometryRepresentation<DEBUG_LEVEL, NDIM>*>> m_bodyMap;
  std::vector<MinMaxType>                                                              m_elementMinMax;
  KDTree<DEBUG_LEVEL, NDIM>                                                            m_kd;
};

#endif // GRIDGENERATOR_GEOMETRY_H
