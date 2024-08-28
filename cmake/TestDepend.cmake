# 引入doctest
include(FetchContent)
FetchContent_Declare(
        doctest
        GIT_REPOSITORY https://gitee.com/acking-you/doctest.git
        GIT_TAG v2.4.9
        GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(doctest)


#FetchContent_Declare(
#  nanobench
#  GIT_REPOSITORY https://gitee.com/tjopenlab/nanobench.git
#  GIT_TAG v4.3.11
#  GIT_SHALLOW TRUE)
#
#FetchContent_MakeAvailable(nanobench)
