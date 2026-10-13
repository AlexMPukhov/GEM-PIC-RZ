# ctest driver for tests/regression.sh (called with -DEXE, -DCMP, -DSRC)
# optional -DBACKEND=serial|openmp|cuda
if(DEFINED BACKEND)
  set(extra numerics.backend=${BACKEND})
endif()
execute_process(COMMAND sh ${SRC}/tests/regression.sh ${EXE} ${CMP} ${extra} RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "regression test failed")
endif()
