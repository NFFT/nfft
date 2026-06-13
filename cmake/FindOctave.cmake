# FindOctave — locate Octave for building .mex modules. Sets:
#   Octave_FOUND, Octave_CONFIG/Octave_MKOCTFILE/Octave_CLI,
#   Octave_INCLUDE_DIRS, Octave_LIBRARY_DIRS, Octave_LIBRARIES,
#   Octave_MEX_EXTENSION ("mex"), Octave_LIBOCTMEX_SOVERSION ("" if none).
find_program(Octave_CONFIG    NAMES octave-config HINTS ${OCTAVE_ROOT} PATH_SUFFIXES bin)
find_program(Octave_MKOCTFILE NAMES mkoctfile     HINTS ${OCTAVE_ROOT} PATH_SUFFIXES bin)
find_program(Octave_CLI       NAMES octave-cli octave HINTS ${OCTAVE_ROOT} PATH_SUFFIXES bin)

set(Octave_MEX_EXTENSION "mex")
set(Octave_LIBOCTMEX_SOVERSION "")

if(Octave_CONFIG)
  execute_process(COMMAND ${Octave_CONFIG} --print OCTINCLUDEDIR
    OUTPUT_VARIABLE Octave_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${Octave_CONFIG} --print OCTLIBDIR
    OUTPUT_VARIABLE Octave_LIBRARY_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

set(Octave_LIBRARIES "")
if(Octave_MKOCTFILE)
  execute_process(COMMAND ${Octave_MKOCTFILE} -p LIBOCTMEX
    OUTPUT_VARIABLE _liboctmex OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${Octave_MKOCTFILE} -p LIBOCTINTERP
    OUTPUT_VARIABLE _liboctinterp OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${Octave_MKOCTFILE} -p LIBOCTAVE
    OUTPUT_VARIABLE _liboctave OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_liboctmex)
    set(Octave_LIBRARIES "${_liboctmex}")
    file(GLOB _octmex_so "${Octave_LIBRARY_DIRS}/liboctmex.so.*"
                         "${Octave_LIBRARY_DIRS}/liboctmex.*.dylib")
    if(_octmex_so)
      list(GET _octmex_so 0 _first)
      string(REGEX MATCH "liboctmex\\.so\\.([0-9]+)" _m "${_first}")
      if(CMAKE_MATCH_1)
        set(Octave_LIBOCTMEX_SOVERSION "${CMAKE_MATCH_1}")
      endif()
    endif()
  else()
    set(Octave_LIBRARIES "${_liboctinterp} ${_liboctave}")  # older Octave: mex in liboctinterp
  endif()
  separate_arguments(Octave_LIBRARIES UNIX_COMMAND "${Octave_LIBRARIES}")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Octave
  REQUIRED_VARS Octave_MKOCTFILE Octave_INCLUDE_DIRS Octave_LIBRARY_DIRS Octave_CLI)
