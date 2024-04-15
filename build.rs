fn main() {
    cxx_build::bridge("src/lib.rs")
        .file("src/geogram_ffi.cpp")
        .file("include/geogram_predicates_psm/Predicates_psm.cpp") // we need to add the ..._psm.cpp to the compile list, just as when compiling in c++
        .std("c++20")
        .compile("cxx-lab");

    println!("cargo:rerun-if-changed=src/lib.rs");
    println!("cargo:rerun-if-changed=src/geogram_ffi.cpp");
    println!("cargo:rerun-if-changed=include/geogram_ffi.h");
}
