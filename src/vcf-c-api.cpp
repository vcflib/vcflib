/*
  C API provides an application binary interface for external use
 */

extern "C" {
#include "vcf-c-api.h"
}

#include "Variant.h"

using namespace std;
using namespace vcflib;

void testme() {
}

void *zig_variant(void *var) {
    return 0L;
}

void *var_parse(const char *line, bool parse_samples) {
    cerr << "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" << endl;
    Variant * var = new Variant(); // construct buffer
    // Variant::parse(string& line, bool parseSamples) {
    string s = line;
    var->parse(s, parse_samples);
    cerr << "HEY\n" << s << "{" << var->id << "}" << endl;
    printf("<%p %s>\n",var,var->id.c_str());
    return var;
}

const char *var_id(void *var) {
    auto v = static_cast<Variant*>(var);
    // cout << "BACK IN C++ getname " << v->id << endl;
    return (v->id.data());
}

const long var_pos(void *var) {
    return (static_cast<Variant*>(var)->position);
}

const char *var_ref(void *var) {
    auto v = static_cast<Variant*>(var);
    return (v->ref.data());
}

const unsigned long var_alt_num(void *variant) {
    auto v = static_cast<Variant*>(variant);
    return v->alt.size();
}

const char **var_alt(void *variant,const char **ret) {
    auto v = static_cast<Variant*>(variant);
    int idx = 0;
    for (auto &a: v->alt) {
    // auto a = &v->alt[0];
        ret[idx] = a.c_str();
        // printf("ptr=%p:%s,%p:%s:\n",a.c_str(),a.c_str(),ret[idx],ret[idx]);
        // printf("ret=%s,%p,%s\n",ret,*ret,*ret);
        // break;
        idx++;
    }
    return ret;
}

void var_set_id(void *var, const char *id) {
    auto v = static_cast<Variant*>(var);
    v->id = id; // copies content
}

void var_set_ref(void *var, const char *ref) {
    auto v = static_cast<Variant*>(var);
    v->ref = ref; // copies content
}

void var_clear_alt(void *var) {
    auto v = static_cast<Variant*>(var);
    v->alt.clear();
}

void var_set_alt(void *var, const char *alt, long idx) {
    auto v = static_cast<Variant*>(var);
    // v->ref = ref; // copies content
    // printf("[C] %s\n",alt);
    v->alt.push_back(alt);
}
