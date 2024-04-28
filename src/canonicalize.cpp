
#include "canonicalize.h"

namespace vcflib {

bool VariantCanonical::canonicalize(FastaReference& fasta_reference, vector<FastaReference*> insertions, bool place_seq, int min_size){

    // Nobody should call this without checking
    assert(canonicalizable());

    // Nobody should call this twice
    assert(!this->canonical);

    // Find where the inserted sequence can come from for insertions
    bool do_external_insertions = !insertions.empty();
    FastaReference* insertion_fasta;
    if (do_external_insertions){
        insertion_fasta = insertions[0];
    }

    bool ref_valid = allATGCN(ref);

    if (!ref_valid && !place_seq){
        // If the reference is invalid, and we aren't allowed to change the ref sequence,
        // we can't canonicalize the variant.
        return false;
    }

    // Check the alts to see if they are not symbolic
    vector<bool> alt_i_atgcn (alt.size());
    for (int i = 0; i < alt.size(); ++i){
        alt_i_atgcn[i] = allATGCN(alt[i]);
    }

    // Only allow single-alt variants
    bool single_alt = alt.size() == 1;
    if (!single_alt){
        // TODO: this will need to be remove before supporting multiple alleles
        cerr << "Warning: multiple ALT alleles not yet allowed for SVs" << endl;
        return false;
    }

    // Fill in the SV tags
    string svtype = getSVTYPE();
    bool has_len = this->info.count("SVLEN") && !this->info.at("SVLEN").empty();
    bool has_seq = this->info.count("SEQ") && !this->info.at("SEQ").empty();
    bool has_span = this->info.count("SPAN") && !this->info.at("SPAN").empty();
    bool has_end = this->info.count("END") && !this->info.at("END").empty();

    // Where is the end, or where should it be?
    long info_end = 0;
    if (has_end) {
        // Get the END from the tag
        info_end = stol(this->info.at("END")[0]);
    }
    else if(ref_valid && !place_seq) {
        // Get the END from the reference sequence, which is ready.
        info_end = this->position + this->ref.length() - 1;
    }
    else if ((svtype == "DEL" || svtype == "INV") && has_span) {
        // For deletions and inversions, we can get the END from the SPAN
        info_end = this->position + abs(stol(this->info.at("SPAN")[0]));
    }
    else if (svtype == "DEL" && has_len) {
        // For deletions, we can get the END from the SVLEN
        info_end = this->position + abs(stol(this->info.at("SVLEN")[0]));
    }
    else if (svtype == "INS"){
        // For insertions, END is just POS if not specified
        info_end = this->position;
    }
    else{
        cerr << "Warning: could not set END info " << *this << endl;
        return false;
    }

    // Commit back the END
    this->info["END"].resize(1);
    this->info["END"][0] = to_string(info_end);
    has_end = true;

    // What is the variant length change?
    // We store it as absolute value
    long info_len = 0;
    if (has_len){
        // Get the SVLEN from the tag
        info_len = abs(stol(this->info.at("SVLEN")[0]));
    }
    else if ((svtype == "INS" || svtype == "DEL") && has_span){
        info_len = abs(stol(this->info.at("SPAN")[0]));
    }
    else if (svtype == "DEL"){
        // We always have the end by now
        // Deletion ends give you length change
        info_len = info_end - this->position;
    }
    else if (svtype == "INV"){
        // Inversions have 0 length change unless otherwise specified.
        info_len = 0;
    }
    else if (svtype == "INS" && has_seq) {
        // Insertions can let us pick it up from the SEQ tag
        info_len = this->info.at("SEQ").at(0).size();
    }
    else{
        cerr << "Warning: could not set SVLEN info " << *this << endl;
        return false;
    }

    // Commit the SVLEN back
    if (svtype == "DEL"){
        // Should be saved as negative
        this->info["SVLEN"].resize(1);
        this->info["SVLEN"][0] = to_string(-info_len);
    }
    else{
        // Should be saved as positive
        this->info["SVLEN"].resize(1);
        this->info["SVLEN"][0] = to_string(info_len);
    }
    // Now the length change is known
    has_len = true;

    // We also compute a span
    long info_span = 0;
    if (has_span){
        // Use the specified span
        info_span = abs(stol(this->info.at("SVLEN")[0]));
    }
    else if (svtype == "INS" || svtype == "DEL"){
        // has_len is always true here
        // Insertions and deletions let us determine the span from the length change, unless they are complex.
        info_span = info_len;
    }
    else if (svtype == "INV"){
        // has_end is always true here
        // Inversion span is start to end
        info_span = info_end - this->position;
    }
    else{
        cerr << "Warning: could not set SPAN info " << *this << endl;
        return false;
    }

    // Commit the SPAN back
    this->info["SPAN"].resize(1);
    this->info["SPAN"][0] = to_string(info_span);
    // Now the span change is known
    has_span = true;

    if (info_end < this->position) {
        cerr << "Warning: SV END is before POS [canonicalize] " <<
        *this << endl << "END: " << info_end << "  " << "POS: " << this->position << endl;
        return false;
    }

    if (has_seq) {
        // Force the SEQ to upper case, if already present
        this->info["SEQ"].resize(1);
        this->info["SEQ"][0] = toUpper(this->info["SEQ"][0]);
    }

    // Set the other necessary SV Tags (SVTYPE, SEQ (if insertion))
    // Also check for agreement in the position tags
    if (svtype == "INS"){
        if (info_end != this->position){
            cerr << "Warning: insertion END and POS do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
            *this << endl << "END: " << info_end << "  " << "POS: " << this->position << endl;

            if (info_end == this->position + info_len) {
                // We can probably guess what they meant here.
                cerr << "Warning: VCF writer incorrecty produced END = POS + SVLEN for an insertion. Fixing END to POS." << endl;
                info_end = this->position;
                this->info["END"][0] = to_string(info_end);
            } else {
                return false;
            }
        }

        if (info_len != info_span){
            cerr << "Warning: insertion SVLEN and SPAN do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
            *this << endl << "SVLEN: " << info_len << "  " << "SPAN: " << info_span << endl;
            return false;
        }

        if (has_seq && allATGCN(this->info.at("SEQ")[0]) && this->info.at("SEQ")[0].size() != info_len){
            cerr << "Warning: insertion SVLEN and SEQ do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
            *this << endl << "SVLEN: " << info_len << "  " << "SEQ length: " << this->info.at("SEQ")[0].size() << endl;
            return false;
        }

        // Set REF
        string ref_base = toUpper(fasta_reference.getSubSequence(this->sequenceName, this->zeroBasedPosition(), 1));
        if (place_seq){
            this->ref.assign(ref_base);
        }

        if (has_seq &&
                 alt[0] != this->info.at("SEQ")[0] &&
                 allATGCN(this->info.at("SEQ")[0])){
            // Try to remove prepended ref sequence, assuming it's left-aligned
            string s = this->alt[0];
            s = toUpper(s.substr(this->ref.length()));
            if (s != this->info.at("SEQ")[0] && !place_seq){
                cerr << "Warning: INS sequence in alt field does not match SEQ tag" << endl <<
                this->alt[0] << " " << this->info.at("SEQ")[0] << endl;
                return false;
            }
            if (place_seq){
                this->alt[0].assign( ref_base + this->info.at("SEQ")[0] );
            }

        }
        else if (alt_i_atgcn[0] && !has_seq){
            string s = this->alt[0];
            s = toUpper(s.substr(this->ref.length()));
            this->info["SEQ"].resize(1);
            this->info.at("SEQ")[0].assign(s);

            if (s.size() != info_len){
                cerr << "Warning: insertion SVLEN and added bases do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
                *this << endl << "SVLEN: " << info_len << "  " << "added bases: " << s.size() << endl;
                return false;
            }

        }
        else if (alt[0][0] == '<' && do_external_insertions){

            string ins_seq;
            string seq_id = alt[0].substr(1, alt[0].size() - 2);

            if (insertion_fasta->index->find(seq_id) != insertion_fasta->index->end()){
                ins_seq = toUpper(insertion_fasta->getSequence(seq_id));
                if (allATGCN(ins_seq)){
                    this->info["SEQ"].resize(1);
                    this->info["SEQ"][0].assign(ins_seq);
                    if (place_seq){
                        this->alt[0].assign(ref_base + ins_seq);
                    }
                }
                else {
                    cerr << "Warning: Loaded invalid alt sequence for: " << *this << endl;
                    return false;
                }

                if (ins_seq.size() != info_len){
                    cerr << "Warning: insertion SVLEN and FASTA do not agree (complex insertions not canonicalizeable) [canonicalize] " <<
                    *this << endl << "SVLEN: " << info_len << "  " << "FASTA bases: " << ins_seq.size() << endl;
                    return false;
                }
            }
            else{
                cerr << "Warning: could not locate alt sequence for: " << *this << endl;
                return false;
            }

        }
        else{
            cerr << "Warning: could not set SEQ [canonicalize]. " << *this << endl;
            return false;
        }
    }
    else if (svtype == "DEL"){
        // Note that info_len has been abs'd and is always positive
        if (this->position + info_len != info_end){
            cerr << "Warning: deletion END and SVLEN do not agree [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "SVLEN: " << info_len << endl;
            return false;
        }

        if (this->position + info_span != info_end){
            cerr << "Warning: deletion END and SPAN do not agree [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "SPAN: " << info_span << endl;
            return false;
        }

        if (info_end > fasta_reference.sequenceLength(this->sequenceName)) {
            cerr << "Warning: deletion END is past end of sequence [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "length: " << fasta_reference.sequenceLength(this->sequenceName) << endl;
            return false;
        }

        // Set REF
        if (place_seq){
            string del_seq = toUpper(fasta_reference.getSubSequence(this->sequenceName, this->zeroBasedPosition(), info_len + 1));
            string ref_base = toUpper(fasta_reference.getSubSequence(this->sequenceName, this->zeroBasedPosition(), 1));
            this->ref.assign( del_seq );
            this->alt[0].assign( ref_base );
        }
    }
    else if (svtype == "INV"){
        if (this->position + info_span != info_end){
            cerr << "Warning: inversion END and SPAN do not agree [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "SPAN: " << info_span << endl;
            return false;
        }

        if (info_len != 0){
            cerr << "Warning: inversion SVLEN specifies nonzero length change (complex inversions not canonicalizeable) [canonicalize] " <<
            *this << endl << "SVLEN: " << info_len << endl;

            if (info_end == this->position + info_len) {
                // We can probably guess what they meant here.
                cerr << "Warning: VCF writer incorrecty produced END = POS + SVLEN for an inversion. Fixing SVLEN to 0." << endl;
                info_len = 0;
                this->info["SVLEN"][0] = to_string(info_len);
            } else {
                return false;
            }
        }

        if (info_end > fasta_reference.sequenceLength(this->sequenceName)) {
            cerr << "Warning: inversion END is past end of sequence [canonicalize] " << *this << endl <<
            "END: " << info_end << "  " << "length: " << fasta_reference.sequenceLength(this->sequenceName) << endl;
            return false;
        }

        if (place_seq){
            string ref_seq = toUpper(fasta_reference.getSubSequence(this->sequenceName, this->zeroBasedPosition(), info_span + 1));
            // Note that inversions still need an anchoring left base at POS
            string inv_seq = ref_seq.substr(0, 1) + reverse_complement(ref_seq.substr(1));
            this->ref.assign(ref_seq);
            this->alt[0].assign(inv_seq);
        }

    }
    else{
        cerr << "Warning: invalid SV type [canonicalize]:" << *this << endl;
        return false;
    }


    this->updateAlleleIndexes();

    // Check for harmony between ref / alt / tags
    if (this->position > stol(this->info.at("END").at(0))){
        cerr << "Warning: position > END. Possible reference genome mismatch." << endl;
        return false;
    }

    if (svtype == "INS"){
        assert(!this->info.at("SEQ")[0].empty());
    }

    this->canonical = true;
    return true;
}



// To canonicalize a variant, we need either both REF and ALT seqs filled in
// or SVTYPE and SVLEN or END or SPAN or SEQ sufficient to define the variant.
bool VariantCanonical::canonicalizable(){
    bool pre_canon = allATGCN(this->ref);

    for (auto& a : this->alt){
        if (!allATGCN(a)){
            pre_canon = false;
        }
    }

    if (pre_canon){
        // It came in in a fully specified way.
        // TODO: ideally, we'd check to make sure ref/alt lengths, svtypes, and ends line up right here.
        return true;
    }

    string svtype = getSVTYPE();

    if (svtype.empty()){
        // We have no SV type, so we can't interpret things.
        return false;
    }

    // Check the tags
    bool has_len = this->info.count("SVLEN") && !this->info.at("SVLEN").empty();
    bool has_seq = this->info.count("SEQ") && !this->info.at("SEQ").empty();
    bool has_span = this->info.count("SPAN") && !this->info.at("SPAN").empty();
    bool has_end = this->info.count("END") && !this->info.at("END").empty();


    if (svtype == "INS"){
        // Insertions need a SEQ, SVLEN, or SPAN
        return has_seq || has_len || has_span;
    }
    else if (svtype == "DEL"){
        // Deletions need an SVLEN, SPAN, or END
        return has_len || has_span || has_end;
    }
    else if (svtype == "INV"){
        // Inversions need a SPAN or END
        return has_span || has_end;
    }
    else{
        // Other SV types are unsupported
        // TODO: DUP
        return false;
    }
}


int VariantCanonical::getMaxReferencePos(){
    if (this->canonical && this->info.find("END") != this->info.end()) {
        // We are cannonicalized and must have a correct END

        int end = 0;
        for (auto s : this->info.at("END")){
            // Get the latest one defined.
            end = max(abs(stoi(s)), end);
        }
        // Convert to 0-based.
        return end - 1;

    }

    if (!this->isSymbolicSV()){
        // We don't necessarily have an END, but we don't need one
        return this->zeroBasedPosition() + this->ref.length() - 1;
    }

    if (this->canonicalizable()){
        // We aren't canonical, but we could be.
        if (this->info.find("END") != this->info.end()){
            // We have an END; blindly trust it
            int end = 0;
            for (auto s : this->info.at("END")){
                // Get the latest one defined.
                end = max(abs(stoi(s)), end);
            }
            // Convert to 0-based.
            return end - 1;

        }
        else if (this->info.find("SVLEN") != this->info.end()){
            // There's no endpoint, but we know an SVLEN.
            // A negative SVLEN means a deletion, so if we find one we can say we delete that much.
            int deleted = 0;
            for (auto s : this->info.at("SVLEN")){
                int alt_len = stoi(s);
                if (alt_len > 0){
                    // Not a deletion, so doesn't affect any ref bases
                    continue;
                }
                deleted = max(-alt_len, deleted);
            }

            // The anchoring base at POS gets added in (because it isn't
            // deleted) but then subtracted out (because we have to do that to
            // match non-SV deletions). For insertions, deleted is 0 and we
            // return 0-based POS. Inversions must have an END.
            return this->zeroBasedPosition() + deleted;
        }
        else{
            cerr << "Warning: insufficient length information for " << *this << endl;
            return -1;
        }
    }
    else {
        cerr << "Warning: can't get end of non-canonicalizeable variant " << *this << endl;
    }
    return -1;
}

} // namespace vcflib
