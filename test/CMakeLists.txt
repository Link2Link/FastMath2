
project(project_test)
message(STATUS "---------------------------")
message(STATUS "Current : ${PROJECT_NAME}")



SET(TEST_INSTALL_DIR ${PROJECT_TOP_DIR}/bin/test)
include(${PROJECT_TOP_DIR}/cmake/TestDepend.cmake)

enable_testing()

FILE(GLOB TEST_SRC "${PROJECT_SOURCE_DIR}/*.cpp")

foreach(item IN LISTS TEST_SRC)
    get_filename_component(SRCNAME ${item} NAME_WE) #得到文件的名字，不带扩展名
    add_executable(${SRCNAME} ${item} ${SRC})
    target_link_libraries(${SRCNAME} doctest_with_main)
    target_include_directories(${SRCNAME} PUBLIC ${PROJECT_TOP_DIR}/FastMath)
    add_test(
            NAME ${SRCNAME}
            COMMAND $<TARGET_FILE:${SRCNAME}>
    )
    install(TARGETS ${SRCNAME} DESTINATION ${TEST_INSTALL_DIR})
    message(STATUS "Test : ${SRCNAME} <- ${item}")

    SET(SRCNAMELIST "${SRCNAMELIST} ${SRCNAME}")
endforeach(item)
