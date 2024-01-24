vcpkg_from_github(
    OUT_SOURCE_PATH SOURCE_PATH
    REPO mfornace/openmm-1
    REF 1188a56cce7bfafc2dc8f8e0e6c0d66519dc09fe
    SHA512 3d0b35050ead66113673875cd8d236a948ff2cb0c0ee2aec58fbda3742ee15e56cdb9be79308ac125a206dd722c83ae8de4a1caa98639558af51a7d165d3d30b
    HEAD_REF master
)

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
    OPTIONS
        -DBUILD_TESTING=OFF
        -DCMAKE_CXX_FLAGS="-Ofast -march=native"
        -DCMAKE_BUILD_TYPE=Release
)

vcpkg_install_cmake()

file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/docs")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/docs")

# Handle copyright
file(INSTALL ${SOURCE_PATH}/docs-source/licenses/Licenses.txt DESTINATION ${CURRENT_PACKAGES_DIR}/share/${PORT} RENAME copyright)
