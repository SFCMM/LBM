#ifndef GRIDGENERATOR_IO_H
#define GRIDGENERATOR_IO_H
#include <config.h>
#include <fstream>
#include <iostream>
#include <sfcmm_common.h>
#include "../cell_filter.h"
//#include <csv/csv.hpp>

namespace hidden::_detail {
const std::function<GBool(GInt)> defaultTrue  = [](GInt /*ignored*/) { return true; };
const std::function<GBool(GInt)> defaultFalse = [](GInt /*ignored*/) { return false; };
} // namespace hidden::_detail

// todo: improve
struct IOIndex {
  GString name;
  GString type;
};

/// Namespace containing definitions to write ASCII
namespace ASCII {
// using namespace csv;
using namespace std;

/// Write out a point-based CSV-file.
/// \tparam DIM Dimensionality of the data and points
/// \param fileName Name of the file to write.
/// \param noValues The number of values to write.
/// \param coordinates The coordinates of the points to write.
/// \param index Index for the name of the columns in the csv file.
/// \param values Values to be written to the csv-file.
/// \param cellFilter Filter function to reduce the output points, e.g., only write points on a surface and so on.
template <GInt DIM>
inline void writePointsCSV(const GString& fileName, const GInt noValues, const std::vector<VectorD<DIM>>& coordinates,
                           const CellFilterManager<DIM>* cellFilter, const std::vector<IOIndex>& index = {},
                           const std::vector<std::vector<GString>>& values = {}) {
  ASSERT(index.size() == values.size(), "Invalid values/index size!");

  cerr0 << SP1 << "Writing " << fileName << ".csv" << std::endl;
  logger << SP1 << "Writing " << fileName << ".csv" << std::endl;

  ofstream                      pointFile;
  static constexpr unsigned int N           = 64;
  static constexpr unsigned int buffer_size = 1024 * N;
  std::array<char, buffer_size> buffer{};
  pointFile.rdbuf()->pubsetbuf(buffer.data(), buffer_size);
  pointFile << setprecision(std::numeric_limits<double>::digits10);
  pointFile.open(fileName + ".csv");

  for(GInt id = 0; id < DIM; ++id) {
    pointFile << coordinate::name.at(id);
    if(id + 1 < DIM) {
      pointFile << ",";
    }
  }
  for(const auto& columnHeader : index) {
    pointFile << "," << columnHeader.name;
  }
  pointFile << "\n";

  for(GInt id = 0; id < noValues; ++id) {
    if(!cellFilter->eval(id)) {
      continue;
    }
    const auto& coord = coordinates[id];
    for(GInt i = 0; i < DIM; ++i) {
      pointFile << coord[i];
      if(i + 1 < DIM) {
        pointFile << ",";
      }
    }
    for(const auto& column : values) {
      pointFile << "," << column[id];
    }
    pointFile << "\n";
  }
  pointFile.close();
}

/// Write out a point-based CSV-file.
/// \tparam DIM Dimensionality of the data and points
/// \param fileName Name of the file to write.
/// \param noValues The number of values to write.
/// \param coordinates The coordinates of the points to write.
/// \param index Index for the name of the columns in the csv file.
/// \param values Values to be written to the csv-file.
template <GInt DIM>
inline void writePointsCSV(const GString& fileName, const GInt noValues, const std::vector<VectorD<DIM>>& coordinates,
                           const std::vector<IOIndex>& index = {}, const std::vector<std::vector<GString>>& values = {}) {
  CellFilterManager<DIM> filterAll = CellFilterManager<DIM>();
  writePointsCSV<DIM>(fileName, noValues, coordinates, &filterAll, index, values);
}

} // namespace ASCII

/// Namespace to write ParaView files in  VTK format.
namespace VTK {
using namespace std;

/// Write the VTK header.
/// \return  VTK-header string.
static constexpr auto header() -> string_view {
  return "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
         "  <PolyData>\n";
}

/// Write the VTK footer.
/// \return VTK-footer string.
static constexpr auto footer() -> string_view {
  return "    </Piece>\n"
         "  </PolyData>\n"
         "</VTKFile> \n";
}

/// Header for a piece of the polydata field.
/// \param noOfPoints Number of points to be written in this piece.
/// \return Piece-header string.
static inline auto piece_header(const GInt noOfPoints) -> GString {
  return "<Piece NumberOfPoints=\"" + to_string(noOfPoints)
         + "\" NumberOfVerts=\"1\" NumberOfLines=\"0\" NumberOfStrips=\"0\" "
           "NumberOfPolys=\"0\" > \n";
}

/// Point-field header.
/// \tparam DIM Number of dimensions of the points
/// \tparam BINARY Is this field in binary format?
/// \return Return point-header string.
template <GInt DIM, GBool BINARY = false>
static inline auto point_header() -> GString {
  // todo: it seems that number of components is ignored by paraview
  // todo: allow change of float type
  if(BINARY) {
    //    return "<Points>\n<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"" + to_string(DIM) + "\" format=\"binary\">
    //    \n";
    return "<Points>\n<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"binary\"> \n";
  }
  //  return "<Points>\n<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"" + to_string(DIM) + "\" format=\"ascii\"> \n";
  return "<Points>\n<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\"> \n";
}

/// Footer of the point-field
/// \return String point-field footer.
static constexpr auto point_footer() -> string_view {
  return "        </DataArray>\n"
         "      </Points>\n";
}

/// Vertice field header
/// \tparam BINARY Are the vertices in binary format?
/// \return String of the Vertices-field header.
template <GBool BINARY = false>
static constexpr auto vert_header() -> string_view {
  if(BINARY) {
    return "      <Verts>\n"
           "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"binary\"> \n";
  }
  return "      <Verts>\n"
         "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\"> \n";
}

/// Return the offset-field data header
/// \return String of the offset-field data header.
static constexpr auto offset_data_header() -> string_view {
  return "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\"> \n";
}

/// Footer of the data-header.
/// \return Return the string data-header.
static constexpr auto data_footer() -> string_view { return "        </DataArray> \n"; }

/// Vertice-field footer.
/// \return String of the vertice-field footer.
static constexpr auto vert_footer() -> string_view { return "        </Verts> \n"; }

/// Point-data-header.
/// \return String of the Point-data-header
static constexpr auto point_data_header() -> string_view { return "      <PointData> \n"; }

/// Point-data-array definition for an int32 type.
/// \tparam binary Is this data-array in binary format?
/// \param name Name of the data-array
/// \return String of the data-array header.
template <GBool binary = false>
static inline auto point_data_int32(const GString& name) -> GString {
  if(binary) {
    return "<DataArray type=\"Int32\" Name=\"" + name + "\" format=\"binary\">\n";
  }
  return "<DataArray type=\"Int32\" Name=\"" + name + "\" format=\"ascii\">\n";
}

/// Point-data-array definition for an float64 type.
/// \tparam binary Is this data-array in binary format?
/// \param name Name of the data-array
/// \return String of the data-array header.
template <GBool binary = false>
static inline auto point_data_float64(const GString& name) -> GString {
  if(binary) {
    return "<DataArray type=\"Float64\" Name=\"" + name + "\" format=\"binary\">\n";
  }
  return "<DataArray type=\"Float64\" Name=\"" + name + "\" format=\"ascii\">\n";
}

/// Point data-array footer
/// \return String for the point data-array footer.
static constexpr auto point_data_footer() -> string_view { return "      </PointData> \n"; }

/// Write out a point-based VTK-file in ASCII-format.
/// \tparam DIM Dimensionality of the data and points
/// \tparam BIN Write as a binary file or not
/// \param fileName Name of the file to write.
/// \param maxNoValues The maximum number of values (or cellId) or size of the local memory of coordinates and values.
/// \param coordinates The coordinates of the points to write.
/// \param index Index for the name of the columns in the csv file.
/// \param values Values to be written to the csv-file.
/// \param filter Filter function to reduce the output points, e.g., only write points on a surface and so on.
template <GInt DIM, GBool BIN>
inline void writePoints(const GString& fileName, const GInt maxNoValues, const std::vector<VectorD<DIM>>& coordinates,
                        const CellFilterManager<DIM>* cellFilter, const std::vector<IOIndex>& index = {},
                        const std::vector<std::vector<GString>>& values = {}) {
  TERMM(-1, "Work in progress");
  // number of points to output
  GInt noOutCells = 0;
  for(GInt id = 0; id < maxNoValues; ++id) {
    noOutCells += static_cast<GInt>(cellFilter->eval(id));
  }

  cerr0 << SP1 << "Writing " << fileName << ".vtp with #" << noOutCells << " cells" << std::endl;
  logger << SP1 << "Writing " << fileName << ".vtp with #" << noOutCells << " cells" << std::endl;

  ofstream                                         pointFile;
  static constexpr unsigned int                    no_kbs_buffer = 64;
  static constexpr unsigned int                    buffer_size   = 1024 * no_kbs_buffer;
  static constexpr std::array<std::string_view, 4> padders       = {"", "=", "==", "==="};

  std::array<char, buffer_size> buffer{};
  pointFile.rdbuf()->pubsetbuf(buffer.data(), buffer_size);
  pointFile.open(fileName + ".vtp");

  pointFile << setprecision(std::numeric_limits<double>::digits10);
  pointFile << header();
  pointFile << piece_header(noOutCells);
  pointFile << point_header<DIM, BIN>();

  for(GInt id = 0; id < maxNoValues; ++id) {
    if(!cellFilter->eval(id)) {
      continue;
    }
    const auto& coord = coordinates[id];
    for(GInt i = 0; i < DIM; ++i) {
      pointFile << coord[i];
      if(i + 1 < DIM) {
        pointFile << " ";
      }
    }
    if(DIM < 3) {
      for(GInt i = DIM; i < 3; ++i) {
        pointFile << " 0.0";
      }
    }
    pointFile << "\n";
  }

  pointFile << point_footer();
  pointFile << vert_header();
  for(GInt id = 0; id < noOutCells; ++id) {
    pointFile << id << "\n";
  }
  pointFile << data_footer();
  pointFile << offset_data_header();
  pointFile << noOutCells << "\n";
  pointFile << data_footer();
  pointFile << vert_footer();
  pointFile << point_data_header();
  GInt i = 0;
  for(const auto& column : values) {
    // todo: fix types
    if(index[i].type == "float64" || index[i].type == "float32") {
      pointFile << point_data_float64(index[i].name);
    } else {
      pointFile << point_data_int32(index[i].name);
    }
    i++;
    for(GInt id = 0; id < maxNoValues; ++id) {
      if(!cellFilter->eval(id)) {
        continue;
      }
      pointFile << column[id] << "\n";
    }
    pointFile << data_footer();
  }
  pointFile << point_data_footer();
  pointFile << footer();
}

/// Namespace for functions to write VTK in ASCII format.
namespace ASCII {


// todo: combine to functions ASCII and binary!!

/// Write out a point-based VTK-file in ASCII-format.
/// \tparam DIM Dimensionality of the data and points
/// \param fileName Name of the file to write.
/// \param maxNoValues The maximum number of values (or cellId) or size of the local memory of coordinates and values.
/// \param coordinates The coordinates of the points to write.
/// \param index Index for the name of the columns in the csv file.
/// \param values Values to be written to the csv-file.
/// \param filter Filter function to reduce the output points, e.g., only write points on a surface and so on.
template <GInt DIM>
inline void writePoints(const GString& fileName, const GInt maxNoValues, const std::vector<VectorD<DIM>>& coordinates,
                        const CellFilterManager<DIM>* cellFilter, const std::vector<IOIndex>& index = {},
                        const std::vector<std::vector<GString>>& values = {}) {
  // number of points to output
  GInt noOutCells = 0;
  for(GInt id = 0; id < maxNoValues; ++id) {
    noOutCells += static_cast<GInt>(cellFilter->eval(id));
  }

  cerr0 << SP1 << "Writing " << fileName << ".vtp with #" << noOutCells << " cells" << std::endl;
  logger << SP1 << "Writing " << fileName << ".vtp with #" << noOutCells << " cells" << std::endl;

  ofstream                      pointFile;
  static constexpr unsigned int no_kbs_buffer = 64;
  static constexpr unsigned int buffer_size   = 1024 * no_kbs_buffer;
  std::array<char, buffer_size> buffer{};
  pointFile.rdbuf()->pubsetbuf(buffer.data(), buffer_size);
  pointFile.open(fileName + ".vtp");

  pointFile << setprecision(std::numeric_limits<double>::digits10);
  pointFile << header();
  pointFile << piece_header(noOutCells);
  pointFile << point_header<DIM>();

  for(GInt id = 0; id < maxNoValues; ++id) {
    if(!cellFilter->eval(id)) {
      continue;
    }
    const auto& coord = coordinates[id];
    for(GInt i = 0; i < DIM; ++i) {
      pointFile << coord[i];
      if(i + 1 < DIM) {
        pointFile << " ";
      }
    }
    if(DIM < 3) {
      for(GInt i = DIM; i < 3; ++i) {
        pointFile << " 0.0";
      }
    }
    pointFile << "\n";
  }

  pointFile << point_footer();
  pointFile << vert_header();
  for(GInt id = 0; id < noOutCells; ++id) {
    pointFile << id << "\n";
  }
  pointFile << data_footer();
  pointFile << offset_data_header();
  pointFile << noOutCells << "\n";
  pointFile << data_footer();
  pointFile << vert_footer();
  pointFile << point_data_header();
  GInt i = 0;
  for(const auto& column : values) {
    // todo: fix types
    if(index[i].type == "float64" || index[i].type == "float32") {
      pointFile << point_data_float64(index[i].name);
    } else {
      pointFile << point_data_int32(index[i].name);
    }
    i++;
    for(GInt id = 0; id < maxNoValues; ++id) {
      if(!cellFilter->eval(id)) {
        continue;
      }
      pointFile << column[id] << "\n";
    }
    pointFile << data_footer();
  }
  pointFile << point_data_footer();
  pointFile << footer();
}
} // namespace ASCII

/// Namespace for functions to write VTK in binary format.
namespace BINARY {

// template<typename T>
// inline void writeBinary(ofstream& outstream, const T* data, const GInt length){
//   static constexpr std::array<std::string_view, 4> padders       = {"", "=", "==", "==="};
//
//   const GInt header_val_size = static_cast<GInt>(sizeof(T)) * length;
//   const GInt number_bytes    = binary::BYTE_SIZE + header_val_size;
//   const GInt number_chars    = static_cast<GInt>(gcem::ceil(number_bytes * 8.0 / 6.0));
//   const GInt padding         = 4 - (number_chars % 4);
//   outstream << base64::encodeLE<GInt, 1>(&header_val_size);
//   //      pointFile.write(base64::encodeLE<GInt32, 2>(&tmp_val[0], noOutCells).c_str(), number_chars); //todo: some how broken??
//   //      pointFile << padders[padding];
//   outstream << base64::encodeLE<T, 2>(data, length) << padders[padding];
//   outstream << "\n" << data_footer();
// }

// todo: combine to functions ASCII and binary!!

/// Write out a point-based VTK-file in binary-format.
/// \tparam DIM Dimensionality of the data and points
/// \param fileName Name of the file to write.
/// \param maxNoValues The maximum number of values (or cellId) or size of the local memory of coordinates and values.
/// \param coordinates The coordinates of the points to write.
/// \param index Index for the name of the columns in the csv file.
/// \param values Values to be written to the csv-file.
/// \param filter Filter function to reduce the output points, e.g., only write points on a surface and so on.
template <GInt DIM>
inline void writePoints(const GString& fileName, const GInt maxNoValues, const std::vector<VectorD<DIM>>& coordinates,
                        const CellFilterManager<DIM>* filterFunction, const std::vector<IOIndex>& index = {},
                        const std::vector<std::vector<GString>>& values = {}) {
  // number of points to output
  GInt noOutCells = 0;
  for(GInt id = 0; id < maxNoValues; ++id) {
    noOutCells += static_cast<GInt>(filterFunction->eval(id));
  }

  cerr0 << SP1 << "Writing " << fileName << ".vtp with #" << noOutCells << " cells" << std::endl;
  logger << SP1 << "Writing " << fileName << ".vtp with #" << noOutCells << " cells" << std::endl;

  ofstream                                         pointFile;
  static constexpr unsigned int                    no_kbs_buffer = 64;
  static constexpr unsigned int                    buffer_size   = 1024 * no_kbs_buffer;
  static constexpr std::array<std::string_view, 4> padders       = {"", "=", "==", "==="};

  std::array<char, buffer_size> buffer{};
  pointFile.rdbuf()->pubsetbuf(buffer.data(), buffer_size);
  pointFile.open(fileName + ".vtp");

  pointFile << header();
  pointFile << piece_header(noOutCells);
  pointFile << point_header<DIM, true>();

  {
    GInt                actualValues = 0;
    std::vector<GFloat> tmp_coords(noOutCells * 3);

    for(GInt id = 0; id < maxNoValues; ++id) {
      if(!filterFunction->eval(id)) {
        continue;
      }
      const auto& coord = coordinates[id];
      for(GInt i = 0; i < DIM; ++i) {
        tmp_coords[actualValues++] = coord[i];
      }
      // pad to 3D
      for(GInt i = DIM; i < 3; ++i) {
        tmp_coords[actualValues++] = 0.0;
      }
    }

    const GInt header_coord_size = 4 * actualValues;
    const GInt number_bytes      = 8 + header_coord_size;
    const GInt number_chars      = static_cast<GInt>(gcem::ceil(static_cast<GDouble>(number_bytes) * 8.0 / 6.0));
    const GInt padding           = 4 - (number_chars % 4);
    pointFile << base64::encodeLE_header<GFloat, GInt>(tmp_coords.data(), actualValues) << padders[padding];
  }

  pointFile << "\n" << point_footer();
  pointFile << vert_header<true>();
  {
    std::vector<GInt> tmp_id(noOutCells);
    for(GInt id = 0; id < noOutCells; ++id) {
      tmp_id[id] = id;
    }
    const GInt header_coord_size = binary::BYTE_SIZE * noOutCells;
    const GInt number_bytes      = binary::BYTE_SIZE + header_coord_size;
    const GInt number_chars      = static_cast<GInt>(gcem::ceil(number_bytes * 8.0 / 6.0));
    const GInt padding           = 4 - (number_chars % 4);
    pointFile << base64::encodeLE_header<GInt, GInt>(tmp_id.data(), noOutCells) << padders[padding];
  }
  pointFile << "\n";
  pointFile << data_footer();
  pointFile << offset_data_header();
  pointFile << noOutCells << "\n";
  pointFile << data_footer();
  pointFile << vert_footer();
  pointFile << point_data_header();
  {
    GInt i = 0;
    for(const auto& column : values) {
      // todo: fix types
      if(index[i].type == "float64" || index[i].type == "float32") {
        pointFile << point_data_float64<true>(index[i++].name);
      } else {
        pointFile << point_data_int32<true>(index[i++].name);
      }


      //      std::vector<GInt32> tmp_val;
      std::vector<GDouble> tmp_val;
      for(GInt id = 0; id < maxNoValues; ++id) {
        if(filterFunction->eval(id)) {
          //          tmp_val.emplace_back(std::stoi(column[id]));
          tmp_val.emplace_back(std::stod(column[id]));
        }
      }

      //      writeBinary(pointFile, tmp_val.data(), noOutCells);

      const GInt header_val_size = static_cast<GInt>(sizeof(GDouble)) * noOutCells;
      const GInt number_bytes    = binary::BYTE_SIZE + header_val_size;
      const GInt number_chars    = static_cast<GInt>(gcem::ceil(number_bytes * 8.0 / 6.0));
      const GInt padding         = 4 - (number_chars % 4);
      pointFile << base64::encodeLE_header<GDouble, GInt>(tmp_val.data(), noOutCells) << padders[padding];

      pointFile << "\n" << data_footer();
    }
  }
  pointFile << point_data_footer();
  pointFile << footer();
  pointFile.close();
}
} // namespace BINARY
} // namespace VTK

namespace HDF5 {}

#endif // GRIDGENERATOR_IO_H
