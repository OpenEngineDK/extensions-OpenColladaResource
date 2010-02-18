SET(PCRELIB "pcre")
SET(PCRECPPLIB "pcrecpp")
SET(LIBXML2LIB "xml2")

IF(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
  SET(LIBXML2LIB "libxml2_a")
  STRING(COMPARE EQUAL ${CMAKE_BUILD_TYPE} "debug" ISDEBUGENABLED)
  IF(ISDEBUGENABLED)
    SET(PCRELIB "pcre-d")
    SET(PCRECPPLIB "pcrecpp-d")
  ELSE(ISDEBUGENABLED)
    SET(PCRELIB "pcre")
    SET(PCRECPPLIB "pcrecpp")
  ENDIF(ISDEBUGENABLED)
ENDIF(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")

FIND_LIBRARY(LIBPCRE
  NAMES 
  ${PCRELIB}
  PATHS
  ${PROJECT_BINARY_DIR}/lib
  ${PROJECT_SOURCE_DIR}/lib
  ${PROJECT_SOURCE_DIR}/libraries
  ${PROJECT_SOURCE_DIR}/libraries/colladadom2.1/lib
  ${OE_LIB_DIR}/pcre/lib
  ENV LD_LIBRARY_PATH
  ENV LIBRARY_PATH
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  NO_DEFAULT_PATH
)

IF(NOT LIBPCRE)
  MESSAGE ("WARNING: Could not find pcre - depending targets will be disabled.")
ENDIF(NOT LIBPCRE)

FIND_LIBRARY(LIBPCRECPP
  NAMES 
  ${PCRECPPLIB}
  PATHS
  ${PROJECT_BINARY_DIR}/lib
  ${PROJECT_SOURCE_DIR}/lib
  ${PROJECT_SOURCE_DIR}/libraries
  ${PROJECT_SOURCE_DIR}/libraries/colladadom2.1/lib
  ${OE_LIB_DIR}/pcrecpp/lib
  ENV LD_LIBRARY_PATH
  ENV LIBRARY_PATH
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  NO_DEFAULT_PATH
)

IF(NOT LIBPCRECPP)
  MESSAGE ("WARNING: Could not find pcrecpp - depending targets will be disabled.")
ENDIF(NOT LIBPCRECPP)

FIND_LIBRARY(LIBXML2
  NAMES 
  ${LIBXML2LIB}
  PATHS
  ${PROJECT_BINARY_DIR}/lib
  ${PROJECT_SOURCE_DIR}/lib
  ${PROJECT_SOURCE_DIR}/libraries
  ${PROJECT_SOURCE_DIR}/libraries/colladadom2.1/lib
  ENV LD_LIBRARY_PATH
  ENV LIBRARY_PATH
  /usr/lib
  /usr/local/lib
  /opt/local/lib
  NO_DEFAULT_PATH
)

IF(NOT LIBXML2)
  MESSAGE ("WARNING: Could not find libxml2 - depending targets will be disabled.")
ENDIF(NOT LIBXML2)

# Create the extension library
ADD_LIBRARY(Extensions_OpenColladaResource
  Resources/ColladaResource.h
  Resources/ColladaResource.cpp
)

TARGET_LINK_LIBRARIES(Extensions_OpenColladaResource
  OpenEngine_Resources
  # Extension dependencies
  ${OPEN_COLLADA_LIBRARIES}
  ${LIBXML2}
  ${LIBPCRECPP}
  ${LIBPCRE}
)

#MESSAGE(${OPEN_COLLADA_LIBRARIES})
