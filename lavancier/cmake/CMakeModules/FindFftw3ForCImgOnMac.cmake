find_path(FFTW3_INCLUDE_DIR fftw3.h PATHS /opt/local/include /usr/local/include /usr/include)

find_library (FFTW3_LIBRARY
NAMES fftw3
PATHS /usr/lib /usr/local/lib /opt/local/lib)

find_library (FFTW3_THREADS_LIBRARY
NAMES fftw3_threads
PATHS /usr/lib /usr/local/lib /opt/local/lib)

if (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY AND FFTW3_THREADS_LIBRARY)
set (FFTW3_FOUND TRUE)
else (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY AND FFTW3_THREADS_LIBRARY)
set (FFTW3_FOUND FALSE)
endif (FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY AND FFTW3_THREADS_LIBRARY)
if (NOT FFTW3_FOUND)
message (FATAL_ERROR "!!! Fftw3 NOT FOUND !!!")
endif (NOT FFTW3_FOUND)