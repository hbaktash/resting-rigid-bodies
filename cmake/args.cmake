if (TARGET taywee::args)
    return()
endif()

include(CPM)

set(ARGS_BUILD_EXAMPLE CACHE BOOL OFF FORCE)
set(ARGS_BUILD_UNITTESTS CACHE BOOL OFF FORCE)

CPMAddPackage(
    NAME args
    GITHUB_REPOSITORY Taywee/args
    GIT_TAG cc2368ca0d8a962862c96c00fe919e1480050f51
)

set_target_properties(args PROPERTIES FOLDER third_party)
