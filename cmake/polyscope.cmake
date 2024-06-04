if (TARGET polyscope::polyscope)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME polyscope
    GITHUB_REPOSITORY nmwsharp/polyscope
    GIT_TAG v2.2.1
)

set_target_properties(polyscope PROPERTIES FOLDER third_party)
add_library(polyscope::polyscope ALIAS polyscope)
