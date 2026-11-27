# Writes a C++ header that holds one SPIR-V file as a uint32_t array, so the
# solver carries its shaders inside the executable.
# Arguments: -DINPUT=<file.spv> -DOUTPUT=<file.h> -DNAME=<identifier>
file(READ "${INPUT}" _hex HEX)
string(LENGTH "${_hex}" _hex_length)
math(EXPR _words "${_hex_length} / 8")
# A SPIR-V word is stored little endian, so the four bytes are reversed.
string(REGEX REPLACE "(..)(..)(..)(..)" "0x\\4\\3\\2\\1,\n" _body "${_hex}")
file(WRITE "${OUTPUT}"
"// Generated from ${INPUT} by cmake/EmbedSpirv.cmake. Do not edit.
#pragma once
#include <cstddef>
#include <cstdint>
static const uint32_t kSpirv_${NAME}[] = {
${_body}};
static const size_t kSpirvWords_${NAME} = ${_words};
")
