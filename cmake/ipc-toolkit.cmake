# include(FetchContent)
# FetchContent_Declare(
#     ipc_toolkit
#     GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
#     GIT_TAG ${IPC_TOOLKIT_GIT_TAG}
# )
# FetchContent_MakeAvailable(ipc_toolkit)

if (TARGET ipc_toolkit::ipc_toolkit)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    ipc_toolkit
    GIT_REPOSITORY https://github.com/ipc-sim/ipc-toolkit.git
    GIT_TAG main
)
FetchContent_MakeAvailable(ipc_toolkit)