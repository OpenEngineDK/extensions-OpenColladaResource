INCLUDE(${OE_CURRENT_EXTENSION_DIR}/FindOpenCollada.cmake)

IF (OPEN_COLLADA_FOUND) 
  INCLUDE_DIRECTORIES(${OPEN_COLLADA_INCLUDE_DIR})
ELSE (OPEN_COLLADA_FOUND)
  MESSAGE ("WARNING: Could not find OpenCollada - depending targets will be disabled.")
ENDIF (OPEN_COLLADA_FOUND)

