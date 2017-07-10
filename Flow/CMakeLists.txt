# Create a library called "" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} 
${GSL_DIR}/include
${CLHEP_INCLUDE_DIR} 
${CMAKE_SOURCE_DIR}/geobase 
${CMAKE_SOURCE_DIR}/parbase
${CMAKE_SOURCE_DIR}/base/source
${CMAKE_SOURCE_DIR}/base/steer
${CMAKE_SOURCE_DIR}/base/event
${CMAKE_SOURCE_DIR}/mpdfield
${CMAKE_SOURCE_DIR}/mcstack
${CMAKE_SOURCE_DIR}/passive
${CMAKE_SOURCE_DIR}/fairtools
${CMAKE_SOURCE_DIR}/tpc
${CMAKE_SOURCE_DIR}/lhetrack
${CMAKE_SOURCE_DIR}/kalman
${CMAKE_SOURCE_DIR}/physics/femto
${CMAKE_SOURCE_DIR}/mpddata
)

include_directories(${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)
 
link_directories(${LINK_DIRECTORIES})

set(PHYS_FLOW_SRC
restore_dca/restore_dca.c
get_dca/get_dca.cxx
get_centrality/get_multiplicity.cxx
create_reduced_tree/reducedTreeCreator.C
)

# fill list of header files from list of source files
# by exchanging the file extension
CHANGE_FILE_EXTENSION(*.cxx *.h PHYS_FLOW_HEADERS "${PHYS_FLOW_SRC}")

set(PHYS_Flow_LINKDEF LinkDef.h)
set(PHYS_FLOW_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/physFlowDict.cxx) 

ROOT_GENERATE_DICTIONARY("${PHYS_FLOW_HEADERS}" "${PHYS_FLOW_LINKDEF}" "${PHYS_FLOW_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

#SET(PHYS_FEMTO_SRCS ${PHYS_FEMTO_SRCS} ${PHYS_FEMTO_DICTIONARY} dbgstream.cxx)
SET(PHYS_FLOW_SRC ${PHYS_FLOW_SRC} ${PHYS_FLOW_DICTIONARY} )

add_library(MpdFlow SHARED ${PHYS_FLOW_SRC})
target_link_libraries(MpdFlow ${ROOT_LIBRARIES})
set_target_properties(MpdFlow PROPERTIES VERSION 0.0.1 SOVERSION 0 )

################ install ###################
install(TARGETS MpdFlow DESTINATION ${CMAKE_BINARY_DIR}/lib)
