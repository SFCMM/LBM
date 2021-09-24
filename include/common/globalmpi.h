// SPDX-License-Identifier: BSD-3-Clause

#ifndef SFCMM_GLOBALMPI_H
#define SFCMM_GLOBALMPI_H

#include <array>
#include <iostream>
#include <mpi.h>

#include "common/sfcmm_types.h"

namespace MPI {

/// Print all information of given MPI_Info object
inline void printMpiInfo(MPI_Info& mpiInfo) {
  int nkeys = 0;

  MPI_Info_get_nkeys(mpiInfo, &nkeys);
  std::cerr << "MPI Info: nkeys = " << nkeys << std::endl;
  for(int i = 0; i < nkeys; i++) {
    std::array<GChar, MPI_MAX_INFO_KEY> key{};
    std::array<GChar, MPI_MAX_INFO_VAL> value{};
    int                                 valuelen = 0;
    int                                 flag     = 0;

    MPI_Info_get_nthkey(mpiInfo, i, &key[0]);
    MPI_Info_get_valuelen(mpiInfo, &key[0], &valuelen, &flag);
    MPI_Info_get(mpiInfo, &key[0], valuelen + 1, &value[0], &flag);
    std::cerr << "MPI Info: [" << i << "] key = " << key.data() << ", value = " << value.data() << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Accessors and storage for global MPI information

/// Class to store global MPI information and to prevent accidental changes
class Information {
 public:
  void init(const GInt domainId, const GInt noDomains) {
    m_globalDomainId  = domainId;
    m_globalNoDomains = noDomains;

    initMPIInformation();
  }

 private:
  void initMPIInformation() {
    MPI_Info_create(&m_mpiInfo);

    // Set header align size to 10KB for netCDF files. Allows to append header data without the need
    // to move all variable data if the header size is exceeded (which may cause MPI I/O errors).
    // Source: https://trac.mcs.anl.gov/projects/parallel-netcdf/wiki/VariableAlignment
    MPI_Info_set(m_mpiInfo, "nc_header_align_size", "10240");
    // Note: possibility to set variable align size
    /* MPI_Info_set(m_mpiInfo, "nc_var_align_size", "4194304"); */


#ifdef MPI_IO_PRINT_INFO
    // Print MPI information on global rank 0
    if(m_globalDomainId == 0) {
      cerr0 << std::endl << "Global MPI information" << std::endl;
      printMpiInfo(m_mpiInfo);
    }
#endif
  }

  friend auto globalDomainId() -> GInt;
  friend auto globalNoDomains() -> GInt;
  friend auto globalMpiInfo() -> const MPI_Info&;
  friend auto isRoot() -> GBool;
  friend auto isSerial() -> GBool;

  GInt     m_globalDomainId  = 0;
  GInt     m_globalNoDomains = 1;
  MPI_Info m_mpiInfo         = MPI_INFO_NULL;
};
inline Information g_mpiInformation; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables)

/// Return global domain id
inline auto globalDomainId() -> GInt { return g_mpiInformation.m_globalDomainId; }
/// Return global number of domains
inline auto globalNoDomains() -> GInt { return g_mpiInformation.m_globalNoDomains; }
/// Return global MPI information
inline auto globalMpiInfo() -> const MPI_Info& { return g_mpiInformation.m_mpiInfo; }
/// Return is root process
inline auto isRoot() -> GBool { return g_mpiInformation.m_globalDomainId == 0; }
/// Return if we running serial
inline auto isSerial() -> GBool { return g_mpiInformation.m_globalNoDomains == 1; }
} // namespace MPI

#endif // SFCMM_GLOBALMPI_H
