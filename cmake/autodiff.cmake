if (TARGET autodiff::autodiff)
    return()
endif()

include(CPM)

set(AUTODIFF_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(AUTODIFF_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(AUTODIFF_BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(AUTODIFF_BUILD_DOCS OFF CACHE BOOL "" FORCE)

CPMAddPackage(
    NAME autodiff
    GITHUB_REPOSITORY autodiff/autodiff
    GIT_TAG v1.1.2
)

set_target_properties(autodiff PROPERTIES FOLDER third_party)
