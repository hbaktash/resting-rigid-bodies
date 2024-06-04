if (TARGET geometry-central::geometry-central)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME geometry-central
    GITHUB_REPOSITORY nmwsharp/geometry-central
    GIT_TAG 095124d31ace21fcda4c06df4b4cb046de84d30a
)

set_target_properties(geometry-central PROPERTIES FOLDER third_party)
add_library(geometry-central::geometry-central ALIAS geometry-central)
