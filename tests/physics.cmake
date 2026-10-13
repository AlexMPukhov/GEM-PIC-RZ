# optional -DBACKEND=serial|openmp|cuda
if(DEFINED BACKEND)
  set(extra numerics.backend=${BACKEND})
endif()
execute_process(COMMAND sh ${SRC}/tests/physics.sh ${EXE} ${CHK} ${extra} RESULT_VARIABLE rc)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "physics test failed")
endif()
