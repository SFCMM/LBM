// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_LOG_H
#define SFCMM_LOG_H
#include <fstream>
#include <map>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <utility>
#include <vector>
#include "macros.h"
#include "util/sys.h"

// todo: add tests

class Log_buffer : public std::stringbuf {
  friend class Log;

 public:
  Log_buffer() = default;

  Log_buffer(const GInt argc, GChar** argv) : m_argc(argc), m_argv(argv) {}

  /**
   * \brief Sets the minimum buffer length that has to be reached before the buffer is flushed.
   * \params[in] minFlushSize Minimum buffer length.
   */
  void setMinFlushSize(GInt minFlushSize) { m_minFlushSize = minFlushSize; }

  virtual void close() = 0;

 protected:
  /**
   * \brief Parses the string input and returns the string with XML entities escaped
   * \details This method iterates over each character of the given input string str and replaces relevant XML
   *          entities with their escaped counterparts.
   *          This code is adapted from http://www.mdawson.net/misc/xmlescape.php (Matson Dawson, 2009).
   *
   * \param[in] str Input string that has characters which need escaping.
   * \return Modified string with XML entities escaped.
   */
  static auto encodeXml(const GString& inputStr) -> GString {
    std::ostringstream tmpEncodeBuffer; // Used as a temporary string buffer

    // Create a for loop that uses an iterator to traverse the complete string
    for(GString::const_iterator iter = inputStr.begin(); iter < inputStr.end(); iter++) {
      // Get current character
      auto c = static_cast<GChar>(*iter);

      // Use a switch/case statement for the five XML entities
      switch(c) {
        case '"':
          tmpEncodeBuffer << "&quot;";
          break; // Replace double quotes
        case '&':
          tmpEncodeBuffer << "&amp;";
          break; // Replace ampersand
        case '\'':
          tmpEncodeBuffer << "&apos;";
          break; // Replace single quote
        case '<':
          tmpEncodeBuffer << "&lt;";
          break; // Replace less-than sign
        case '>':
          tmpEncodeBuffer << "&gt;";
          break; // Replace greater-than sign
        default:
          tmpEncodeBuffer << c; // By default, just append current character
      }
    }

    // Return encoded stream as a string
    return tmpEncodeBuffer.str();
  }

  /**
   * \brief Creates an XML prefix using the domain id that is prepended to each message.
   * Makes use of an array attribute filled by the user to generate the XML string.
   *
   */
  virtual void createPrefixMessage() {
    // Create temporary stream
    std::ostringstream tmpStream;

    // Fill stream with formatted domain id
    tmpStream << "<m d=\"" << m_domainId << "\" ";

    for(const auto& attribute : m_prefixAttributes) {
      tmpStream << attribute.first << "=\"" << attribute.second() << "\" ";
    }

    tmpStream << ">";

    // Set prefix message to tmpBuffer string
    m_prefixMessage = tmpStream.str();
  }

  /**
   * \brief Creates an XML suffix that is appended to each message.
   */
  virtual void createSuffixMessage() { m_suffixMessage = "</m>\n"; }


  /**
   * \brief Return an XML header that should written at the beginning of each log file.
   * \return The XML header.
   */
  virtual auto getXmlHeader() -> GString {
    using namespace std;

    // Gets the current executionCommand
    stringstream executionCommand;
    executionCommand.str("");
    executionCommand << m_argv[0];
    for(GInt n = 1; n < m_argc; n++) {
      executionCommand << " " << m_argv[n];
    }

    // Create temporary buffer
    ostringstream tmpBuffer;

    // Write XML header information to buffer
    tmpBuffer << R"(<?xml version="1.0" standalone="yes" ?>)"
              << "\n";
    tmpBuffer << R"(<root>)"
              << "\n";
    tmpBuffer << R"(<meta name="noDomains" content=")" << m_noDomains << "\" />\n";
    tmpBuffer << R"(<meta name="dateCreation" content=")" << dateString() << "\" />\n";
    tmpBuffer << R"(<meta name="fileFormatVersion" content=")" << m_fileFormatVersion << "\" />\n";
    tmpBuffer << R"(<meta name="user" content=")" << userString() << "\" />\n";
    tmpBuffer << R"(<meta name="host" content=")" << hostString() << "\" />\n";
    tmpBuffer << R"(<meta name="dir" content=")" << getCWD() << "\" />\n";
    tmpBuffer << R"(<meta name="executionCommand" content=")" << executionCommand.str() << "\" />\n";
    tmpBuffer << R"(<meta name="revision" content=")" << XSTRINGIFY(PROJECT_VER) << "\" />\n";
    tmpBuffer << R"(<meta name="build" content=")" << XSTRINGIFY(COMPILER_NAME) << " " << XSTRINGIFY(BUILD_TYPE) << " ("
              << GString(XSTRINGIFY(COMPILER_VER)) << ")"
              << "\" />\n";

    // Return XML header
    return tmpBuffer.str();
  }

  /**
   * \brief Return an XML footer that should written at the end of each log file.
   * \return The XML footer.
   */
  virtual auto getXmlFooter() -> GString {
    // Create temporary buffer
    std::ostringstream tmpBuffer;

    // Write XML footer to buffer
    tmpBuffer << R"(<meta name="dateClosing" content=")" << dateString() << "\" />\n";
    tmpBuffer << "</root>\n";

    // Return XML footer
    return tmpBuffer.str();
  }

  /// Flush the current buffer to output target
  virtual void flushBuffer() = 0;

  /// Current domain Id
  /// \return domain Id
  [[nodiscard]] auto inline domainId() const -> int { return m_domainId; }

  /// Current domain Id (reference version)
  /// \return reference to domain Id
  auto inline domainId() -> int& { return m_domainId; }

  /// Current number of MPI ranks
  /// \return number of MPI ranks
  [[nodiscard]] auto inline noDomains() const -> int { return m_noDomains; }

  /// Current number of MPI ranks (reference version)
  /// \return reference to the number of MPI ranks
  auto inline noDomains() -> int& { return m_noDomains; }

  /// Prefix message of the XML output of the log message
  /// \return Prefix message of the XML output
  [[nodiscard]] auto inline prefixMessage() const -> const GString& { return m_prefixMessage; }

  /// Suffix message of the XML output of the log message
  /// \return Suffix message of the XML output
  [[nodiscard]] auto inline suffixMessage() const -> const GString& { return m_suffixMessage; }
  [[nodiscard]] auto inline minFlushSize() const -> GInt { return m_minFlushSize; }

 private:
  static constexpr GInt m_fileFormatVersion = 1; //!< File format version (increase this by one every time you make
  //!< changes that could affect postprocessing tools)
  int     m_domainId{0};     //!< Contains the MPI rank (= domain id) of this process
  int     m_noDomains{1};    //!< Contains the MPI rank count (= number of domains)
  GInt    m_minFlushSize{0}; //!< Minimum length of the internal buffer before flushing
  GString m_prefixMessage;   //!< Stores the prefix that is prepended to each output
  GString m_suffixMessage;   //!< Stores the suffix that is appended to each output
  GInt    m_argc{};          /// Command line argument count
  GChar** m_argv{};          /// Command line arguments

  std::map<GString, std::function<GString()>> m_prefixAttributes; /// attributes of the XML log message
};


/**
 * \brief Customized buffer to facilitate of a regular physical file for each processor within an MPI communicator.
 * \details This class can be used as a regular string buffer, as it inherits from stringbuf. On flushing the
 *          buffer, the contents of the buffer are written to a file using an ofstream. This is mainly for cases
 *          where logging speed is crucial, as the implementation is very lightweight and since for each process
 *          an individual file is maintained. There is an option to use this buffer but to create only one file
 *          for the MPI root domain (see #m_rootOnlyHardwired).
 */
class Log_simpleFileBuffer : public Log_buffer {
 public:
  Log_simpleFileBuffer() = default;

  /**
   * \brief This constructor yields a new instance that can immediately be used to write messages to a regular file.
   * \details Internally, this constructor just passes the parameters to open() (see open() for more details).
   *
   * \param[in] filename          Filename that should be used for the file.
   * \param[in] mpiComm           MPI communicator which is used to determine rank/domain information.
   * \param[in] rootOnlyHardwired If true, only rank 0 creates a file and writes to it. On all other processors, no
   *                              file is opened and at each flushing of the buffer, the buffer content is discarded.
   */
  Log_simpleFileBuffer(const GString& filename, const GInt argc, GChar** argv, MPI_Comm mpiComm = MPI_COMM_WORLD,
                       GBool rootOnlyHardwired = false)
    : Log_buffer(argc, argv) {
    open(filename, mpiComm, rootOnlyHardwired);
  }

  ~Log_simpleFileBuffer() override { Log_simpleFileBuffer::close(); };

  Log_simpleFileBuffer(const Log_simpleFileBuffer&) = delete;
  Log_simpleFileBuffer(Log_simpleFileBuffer&&)      = delete;
  auto operator=(const Log_simpleFileBuffer&) -> Log_simpleFileBuffer& = delete;
  auto operator=(Log_simpleFileBuffer&&) -> Log_simpleFileBuffer& = delete;

  /**
   * \brief Initialization of the file I/O environment.
   * \details After a successful call to this method the file stream is ready to use. Any previous information
   *          written to the buffer is lost when open is called. This function creates a new file as specified in
   *          filename (deleting any existing files with the same name), and creates the XML prefixes and suffixes
   *          used for each message. It also writes the necessary XML header information to the file.
   *
   * \param[in] filename          Filename that should be used for the file.
   * \param[in] mpiComm           MPI communicator which is used to determine rank/domain information.
   * \param[in] rootOnlyHardwired If true, only rank 0 creates a file and writes to it. On all other processors, no
   *                              file is opened and at each flushing of the buffer, the buffer content is discarded.
   */
  void open(const GString& filename, MPI_Comm mpiComm = MPI_COMM_WORLD, GBool rootOnlyHardwired = false) {
    // Open file only if it was not yet done
    if(!m_isOpen) {
      // Set MPI communicator group
      m_mpiComm = mpiComm;

      // Get domain id and number of domains
      MPI_Comm_rank(m_mpiComm, &domainId());
      MPI_Comm_size(m_mpiComm, &noDomains());

      // Set whether only domain 0 should do any writing (including the creation of a file)
      m_rootOnlyHardwired = rootOnlyHardwired;

      // Only open the file if m_rootOnlyHardwired was not set as true. Otherwise, the file state remains closed.
      if(!(m_rootOnlyHardwired && domainId() != 0)) {
        // Set filename
        m_filename = filename;

        // Open file
        m_file.open(m_filename.c_str());

        // Clear internal buffer in order to dismiss any previous input
        str("");

        // Create prefix and suffix messages
        createPrefixMessage();
        createSuffixMessage();

        // Write root and meta information to file
        m_file << getXmlHeader() << std::flush;

        // Set state variable
        m_isOpen = true;
      }
    }
  }

  /**
   * \brief Closes the file.
   * \details Any subsequent write statements to the file stream are discarded after this method is called. After
   *          close() is called, an XML footer is written to the file. Then the file is closed.
   */
  void close() override {
    // Only close file if was opened before
    if(m_isOpen) {
      // Force flushing of the internal buffer
      flushBuffer();

      // Write XML footer to file and flush stream
      m_file << getXmlFooter() << std::flush;

      // Close file stream
      m_file.close();

      // Set state variable
      m_isOpen = false;
    }
  }

 protected:
  /**
   * \brief Flushes the buffer by writing the contents to the file.
   * \details Sync is called automatically when an "endl" is sent to the stream. At first the buffer content is
   *          wrapped in the prefix and suffix messages, then the entire string is written to the file by calling
   *          flushBuffer(). Finally, the internal buffers are reset.
   *
   * \return Zero by default.
   */
  auto sync() -> int override {
    // Only write if the file was already opened
    if(m_isOpen) {
      // Create formatted string, escape any XML entities in the message, and save to temporary buffer
      m_tmpBuffer << prefixMessage() << encodeXml(str()) << suffixMessage();

      // Only write to file if current buffer length exceeds the minimum size for flushing
      if(m_tmpBuffer.str().length() >= static_cast<unsigned>(minFlushSize())) {
        // Write the string to the file and flush the stream
        m_file << m_tmpBuffer.str() << std::flush;

        // Reset temporary buffer
        m_tmpBuffer.str("");
      }
    }

    // Reset internal buffer
    str("");

    // Default return value for sync()
    return 0;
  }

  /**
   * \brief Flushes the buffer by writing the contents to the file.
   * \details Sync is called automatically when an "endl" is sent to the stream. At first the buffer content is
   *          wrapped in the prefix and suffix messages, then the entire string is written to the file. Finally,
   *          the internal buffers are reset.
   *
   * \return Zero by default.
   */
  void flushBuffer() override {
    // Only write if the file was already opened
    if(m_isOpen) {
      // Write the string to the file and flush the stream
      m_file << m_tmpBuffer.str() << std::flush;

      // Reset temporary buffer
      m_tmpBuffer.str("");
    }
  }

 private:
  GBool              m_isOpen{false};            //!< Stores whether the file(s) were already opened
  GBool              m_rootOnlyHardwired{false}; //!< If true, only domain 0 opens and uses a file
  GString            m_filename;                 //!< Filename on disk
  std::ofstream      m_file;                     //!< File stream tied to physical file on disk
  MPI_Comm           m_mpiComm{};                //!< MPI communicator group
  std::ostringstream m_tmpBuffer;                //!< Temporary buffer to hold string until flushing
};

/**
 * \brief Base class for all Log<xyz> classes.
 * \details This class is used to hold stream/buffer-independent methods. All Log<xyz>
 * subclasses inherit from this class. The auxiliary classes Log_<xyz> (especially the
 * buffers), however, should NOT have this class as their baseclass.
 */
class Log : public std::ostream {
  friend class Log_buffer;

 public:
#ifdef CLANG_COMPILER
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
#endif
  Log() : std::ostream(m_buffer.get()){};
#ifdef CLANG_COMPILER
#pragma clang diagnostic pop
#endif

  /** \brief Adds an attribute to the prefix of the XML string.
   * \param[in] att The attribute to add, consists of a pair of MStrings.
   * \return The location of the attribute in the vector of pairs.
   */
  auto addAttribute(const std::pair<GString, std::function<GString()>>& att) -> GInt {
    m_buffer->m_prefixAttributes.emplace(att);
    m_buffer->createPrefixMessage();
    return static_cast<GInt>(m_buffer->m_prefixAttributes.size() - 1);
  }

  /** \brief Erases an attribute from the prefix of the XML string.
   *  \param[in] attName The Name of the attribute to delete.
   */
  void eraseAttribute(const GString& attName) {
    m_buffer->m_prefixAttributes.erase(attName);
    m_buffer->createPrefixMessage();
  }

  /** \brief Modifies an attribute of the prefix of the XML string.
   *  \param[in] attName The Name of the attribute to modify.
   *  \param[in] att The new attribute to replace the old one, given by a pair of MStrings.
   */
  void modifyAttribute(const GString& attName, const std::pair<GString, std::function<GString()>>& att) {
    m_buffer->m_prefixAttributes.erase(attName);
    m_buffer->m_prefixAttributes.emplace(att);
    m_buffer->createPrefixMessage();
  }

  /** \brief Modifies an attribute of the prefix of the XML string.
   *  \param[in] attName The Name of the attribute to modify.
   *  \param[in] att The new attribute to replace the old one, given by a pair of MStrings.
   */
  void updateAttributeValue(const GString& attName, std::function<GString()> att) {
    m_buffer->m_prefixAttributes[attName] = std::move(att);
    m_buffer->createPrefixMessage();
  }

  /** \brief  Updates all attributes of the prefix of the XML string..
   */
  void updateAttributes() { m_buffer->createPrefixMessage(); }

  /// Get access to the buffer of the Log using an unique ptr to log_buffer instance.
  /// \return Unique ptr to the log_buffer instance.
  inline auto buffer() -> std::unique_ptr<Log_buffer>& { return m_buffer; }

  /// Get access to the buffer of the Log using an *const* unique ptr to log_buffer instance.
  /// \return *Const* Unique ptr to the log_buffer instance.
  inline auto buffer() const -> const std::unique_ptr<Log_buffer>& { return m_buffer; }

 private:
  std::unique_ptr<Log_buffer> m_buffer;
};

/**
 * \brief Class to create an output stream for a writable file, using a physical file.
 * \details This class can be used to open a file on all processors (alternatively: only on a specified MPI
 *          communicator) and write to it using the normal C++ stream syntax (i.e. just like cout or cerr).
 *
 *          A regular physical files for each processor, and to write
 *          directly to it using a regular ofstream. This mode also allows for the setting that only process 0
 *          within an MPI communicator opens a file to write to.
 */
class LogFile : public Log {
 public:
  LogFile() = default;
  /**
   * \brief Constructor creates LogFile and calls ostream constructor with reference to it.
   * \details When this constructor is used, the stream is immediately ready to use. For information about the
   * parameters, please have a look at Log*Buffer::open.
   *
   * \param[in] filename Name of the file to open.
   * \param[in] mpiComm MPI communicator for which to open the file.
   */
  explicit LogFile(const GString& filename, MPI_Comm mpiComm = MPI_COMM_WORLD, GBool rootOnlyHardwired = false) {
    open(filename, rootOnlyHardwired, 0, nullptr, mpiComm);
  }
  ~LogFile() override { close(); };

  LogFile(const LogFile&) = delete;
  LogFile(LogFile&&)      = delete;
  auto operator=(const LogFile&) -> LogFile& = delete;
  auto operator=(LogFile&&) -> LogFile& = delete;

  /**
   * \brief Opens a file by passing the parameters to Log_<xyz>FileBuffer::open(...).
   * \details The method then creates a new internal buffer and passes along the parameters.
   *
   * \param[in] filename Name of the file to open.
   * \param[in] mpiComm MPI communicator for which to open the file.
   * \param[in] rootOnlyHardwired If true, only rank 0 creates a file and writes to it. On all other processors, no
   *                              file is opened and at each flushing of the buffer, the buffer content is discarded.
   */
  void open(const GString& filename, const GBool rootOnlyHardwired, const GInt argc, GChar** argv, MPI_Comm mpiComm = MPI_COMM_WORLD) {
    // Only open file if it was not yet opened
    if(!m_isOpen) {
      // Open a simple file
      buffer() = std::make_unique<Log_simpleFileBuffer>(filename, argc, argv, mpiComm, rootOnlyHardwired);

      // Associate the stream with the newly created buffer
      rdbuf(buffer().get());

      // Set state variable
      m_isOpen = true;
    }
  }

  /**
   * \brief Pass the close call to the respective internal buffer.
   * \details All attempts to write to the stream after closing it will be discarded.
   */
  void close() {
    // Only close file if was already opened
    if(m_isOpen) {
      buffer()->close();

      // Delete internal buffer to prevent memory leaks
      buffer().reset();

      // Set state variable
      m_isOpen = false;
    }
  }

  /**
   * \brief Sets the minimum buffer length that has to be reached before the buffer is flushed.
   * \details Flushing the buffer means that the contents of the buffer in memory is written to the file. If the
   *          file stream was not opened yet, this method just returns 0 and does nothing else.
   *
   * \params[in] minFlushSize Minimum buffer length.
   */
  void setMinFlushSize(const GInt minFlushSize) {
    if(buffer() != nullptr) {
      buffer()->setMinFlushSize(minFlushSize);
    }
  }

 private:
  bool m_isOpen = false; //!< Stores whether a file was already opened or not
};
inline LogFile logger; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

#endif // SFCMM_LOG_H
