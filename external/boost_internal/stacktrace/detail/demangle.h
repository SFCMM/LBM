#ifndef LBM_DEMANGLE_H
#define LBM_DEMANGLE_H

inline char const * demangle_alloc( char const * name ) BOOST_NOEXCEPT
{
  int status = 0;
  std::size_t size = 0;
  return abi::__cxa_demangle( name, NULL, &size, &status );
}

class scoped_demangled_name
{
private:
  char const * m_p;

public:
  explicit scoped_demangled_name( char const * name ) BOOST_NOEXCEPT :
      m_p( demangle_alloc( name ) )
  {
  }

  ~scoped_demangled_name() BOOST_NOEXCEPT
  {
    std::free( const_cast< char* >( m_p ) );
  }

  char const * get() const BOOST_NOEXCEPT
  {
    return m_p;
  }

  BOOST_DELETED_FUNCTION(scoped_demangled_name( scoped_demangled_name const& ))
  BOOST_DELETED_FUNCTION(scoped_demangled_name& operator= ( scoped_demangled_name const& ))
};

inline std::string demangle( char const * name )
{
  scoped_demangled_name demangled_name( name );
  char const * p = demangled_name.get();
  if( !p )
    p = name;
  return p;
}
#endif // LBM_DEMANGLE_H
