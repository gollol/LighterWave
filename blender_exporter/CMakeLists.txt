find_package(
    Python3 3.11 EXACT
    COMPONENTS Interpreter
)

if (Python3_Interpreter_FOUND)
    set(output ${CMAKE_BINARY_DIR}/lightwave_blender.zip)

    if(WIN32)
        execute_process(
            COMMAND py -3.11 -c "import cython"
            RESULT_VARIABLE EXIT_CODE
            OUTPUT_QUIET
        )

        if (NOT ${EXIT_CODE} EQUAL 0)
            message(
                FATAL_ERROR
                "The \"cython\" Python3 package is not installed. Please install it using the following command: \"py -3.11 -m pip install cython\"."
            )
        endif()

        execute_process(
            COMMAND py -3.11 -c "import setuptools"
            RESULT_VARIABLE EXIT_CODE
            OUTPUT_QUIET
        )

        if (NOT ${EXIT_CODE} EQUAL 0)
            message(
                FATAL_ERROR
                "The \"setuptools\" Python3 package is not installed. Please install it using the following command: \"py -3.11 -m pip install setuptools\"."
            )
        endif()

        add_custom_target(lightwave_blender ALL
            COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/build_lightwave_blender.bat"
            COMMAND ${CMAKE_COMMAND} -E tar cf "${output}" --format=zip -- "${CMAKE_CURRENT_SOURCE_DIR}/lightwave_blender"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        )
    else()
        execute_process(
            COMMAND python3.11 -c "import cython"
            RESULT_VARIABLE EXIT_CODE
            OUTPUT_QUIET
        )

        if (NOT ${EXIT_CODE} EQUAL 0)
            message(
                FATAL_ERROR
                "The \"cython\" Python3 package is not installed. Please install it using the following command: \"python3.11 -m pip install cython\"."
            )
        endif()

        execute_process(
            COMMAND python3.11 -c "import setuptools"
            RESULT_VARIABLE EXIT_CODES
            OUTPUT_QUIET
        )

        if (NOT ${EXIT_CODE} EQUAL 0)
            message(
                FATAL_ERROR
                "The \"setuptools\" Python3 package is not installed. Please install it using the following command: \"python3.11 -m pip install setuptools\"."
            )
        endif()

        add_custom_target(lightwave_blender ALL
            COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/build_lightwave_blender.sh"
            COMMAND ${CMAKE_COMMAND} -E tar cf "${output}" --format=zip -- "${CMAKE_CURRENT_SOURCE_DIR}/lightwave_blender"
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        )
    endif()
    
    install(FILES "${output}" DESTINATION plugins COMPONENT plugins)
else()
    message(WARNING "Python3.11 not found. Skipping lightwave_blender build.")
endif()
