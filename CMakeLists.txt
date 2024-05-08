cmake_minimum_required(VERSION 3.22.1)
project(Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24)

# Postavljam izlazni direktorij za sve izvršne datoteke
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/build)

# Dodajem podmapu src
#add_subdirectory(src)
add_executable(Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24 main.cpp)
# Dodajem podmapu submodules/bioparser
if (NOT TARGET bioparser)
  add_subdirectory(submodules/bioparser EXCLUDE_FROM_ALL)
endif ()

target_link_libraries(Analiza_metagenomskog_uzorka_na_temelju_distribucije_k-mera-23-24 bioparser) #nadodano linkanje
