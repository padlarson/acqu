# simplistic 

Message(STATUS "Looking for a2display...")

set(A2DISPLAY_SEARCH_PATHS
  ${CMAKE_SOURCE_DIR}/../a2display)

FIND_PATH(A2DISPLAY_INCLUDE_DIR a2display.h
  PATHS ${A2DISPLAY_SEARCH_PATHS}/inc
  NO_DEFAULT_PATH
  )

if(NOT A2DISPLAY_INCLUDE_DIR)
  Message(STATUS "Looking for a2display... - a2display.h not found")
  return()
endif()

FIND_LIBRARY(A2DISPLAY_LIBRARIES NAMES a2display
  PATHS ${A2DISPLAY_SEARCH_PATHS}/build
  NO_DEFAULT_PATH
  )

set(A2DISPLAY_FOUND TRUE)
Message(STATUS "Looking for a2display... - Found ${A2DISPLAY_LIBRARIES}")
