#include <rebind/Document.h>
#include <nupack/Version.h>
#include <boost/exception/diagnostic_information.hpp>
#include <boost/core/demangle.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

namespace nupack {

/// Perform functions declared elsewhere one time; the order should reflect dependencies
void render_submodules(rebind::Document &doc) {
#   define NUPACK_TMP(F) void F(rebind::Document &); F(doc);
    NUPACK_TMP(render_constants);
    NUPACK_TMP(render_math);
    NUPACK_TMP(render_model);
    NUPACK_TMP(render_thermo);
    NUPACK_TMP(render_design);
#   undef NUPACK_TMP
}

/******************************************************************************************/

/// Top level module initializer running at library load time
void write_document(rebind::Document &doc) {
    rebind::set_demangler(&boost::core::demangle);
    try {
        render_submodules(doc);
    } catch (...) {
        std::cerr << "C++ import failed:\n" << boost::current_exception_diagnostic_information() << std::endl;
        throw;
    }
}

/******************************************************************************************/

}

namespace rebind {

void init(Document &doc) {
    nupack::write_document(doc);
}

}
