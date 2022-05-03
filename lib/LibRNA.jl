module LibRNA

using ViennaRNA_jll
export ViennaRNA_jll

using CEnum

struct vrna_md_s
    temperature::Cdouble
    betaScale::Cdouble
    pf_smooth::Cint
    dangles::Cint
    special_hp::Cint
    noLP::Cint
    noGU::Cint
    noGUclosure::Cint
    logML::Cint
    circ::Cint
    gquad::Cint
    uniq_ML::Cint
    energy_set::Cint
    backtrack::Cint
    backtrack_type::Cchar
    compute_bpp::Cint
    nonstandards::NTuple{64, Cchar}
    max_bp_span::Cint
    min_loop_size::Cint
    window_size::Cint
    oldAliEn::Cint
    ribo::Cint
    cv_fact::Cdouble
    nc_fact::Cdouble
    sfact::Cdouble
    rtype::NTuple{8, Cint}
    alias::NTuple{21, Cshort}
    pair::NTuple{21, NTuple{21, Cint}}
    pair_dist::NTuple{7, NTuple{7, Cfloat}}
end

const vrna_md_t = vrna_md_s

struct vrna_basepair_s
    i::Cint
    j::Cint
end

const vrna_basepair_t = vrna_basepair_s

struct vrna_elem_prob_s
    i::Cint
    j::Cint
    p::Cfloat
    type::Cint
end

const vrna_plist_t = vrna_elem_prob_s

struct vrna_bp_stack_s
    i::Cuint
    j::Cuint
end

const vrna_bp_stack_t = vrna_bp_stack_s

struct vrna_cpair_s
    i::Cint
    j::Cint
    mfe::Cint
    p::Cfloat
    hue::Cfloat
    sat::Cfloat
    type::Cint
end

const vrna_cpair_t = vrna_cpair_s

struct vrna_sect_s
    i::Cint
    j::Cint
    ml::Cint
end

const vrna_sect_t = vrna_sect_s

struct vrna_color_s
    hue::Cfloat
    sat::Cfloat
    bri::Cfloat
end

const vrna_color_t = vrna_color_s

struct vrna_data_linear_s
    position::Cuint
    value::Cfloat
    color::vrna_color_t
end

const vrna_data_lin_t = vrna_data_linear_s

const FLT_OR_DBL = Cdouble

const PAIR = vrna_basepair_s

const plist = vrna_elem_prob_s

const cpair = vrna_cpair_s

const sect = vrna_sect_s

const bondT = vrna_bp_stack_s

function vrna_md_set_default(md)
    ccall((:vrna_md_set_default, libRNA), Cvoid, (Ptr{vrna_md_t},), md)
end

function vrna_md_update(md)
    ccall((:vrna_md_update, libRNA), Cvoid, (Ptr{vrna_md_t},), md)
end

function vrna_md_copy(md_to, md_from)
    ccall((:vrna_md_copy, libRNA), Ptr{vrna_md_t}, (Ptr{vrna_md_t}, Ptr{vrna_md_t}), md_to, md_from)
end

function vrna_md_option_string(md)
    ccall((:vrna_md_option_string, libRNA), Ptr{Cchar}, (Ptr{vrna_md_t},), md)
end

function vrna_md_set_nonstandards(md, ns_bases)
    ccall((:vrna_md_set_nonstandards, libRNA), Cvoid, (Ptr{vrna_md_t}, Ptr{Cchar}), md, ns_bases)
end

function vrna_md_defaults_reset(md_p)
    ccall((:vrna_md_defaults_reset, libRNA), Cvoid, (Ptr{vrna_md_t},), md_p)
end

function vrna_md_defaults_temperature(T)
    ccall((:vrna_md_defaults_temperature, libRNA), Cvoid, (Cdouble,), T)
end

function vrna_md_defaults_temperature_get()
    ccall((:vrna_md_defaults_temperature_get, libRNA), Cdouble, ())
end

function vrna_md_defaults_betaScale(b)
    ccall((:vrna_md_defaults_betaScale, libRNA), Cvoid, (Cdouble,), b)
end

function vrna_md_defaults_betaScale_get()
    ccall((:vrna_md_defaults_betaScale_get, libRNA), Cdouble, ())
end

function vrna_md_defaults_pf_smooth(s)
    ccall((:vrna_md_defaults_pf_smooth, libRNA), Cvoid, (Cint,), s)
end

function vrna_md_defaults_pf_smooth_get()
    ccall((:vrna_md_defaults_pf_smooth_get, libRNA), Cint, ())
end

function vrna_md_defaults_dangles(d)
    ccall((:vrna_md_defaults_dangles, libRNA), Cvoid, (Cint,), d)
end

function vrna_md_defaults_dangles_get()
    ccall((:vrna_md_defaults_dangles_get, libRNA), Cint, ())
end

function vrna_md_defaults_special_hp(flag)
    ccall((:vrna_md_defaults_special_hp, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_special_hp_get()
    ccall((:vrna_md_defaults_special_hp_get, libRNA), Cint, ())
end

function vrna_md_defaults_noLP(flag)
    ccall((:vrna_md_defaults_noLP, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_noLP_get()
    ccall((:vrna_md_defaults_noLP_get, libRNA), Cint, ())
end

function vrna_md_defaults_noGU(flag)
    ccall((:vrna_md_defaults_noGU, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_noGU_get()
    ccall((:vrna_md_defaults_noGU_get, libRNA), Cint, ())
end

function vrna_md_defaults_noGUclosure(flag)
    ccall((:vrna_md_defaults_noGUclosure, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_noGUclosure_get()
    ccall((:vrna_md_defaults_noGUclosure_get, libRNA), Cint, ())
end

function vrna_md_defaults_logML(flag)
    ccall((:vrna_md_defaults_logML, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_logML_get()
    ccall((:vrna_md_defaults_logML_get, libRNA), Cint, ())
end

function vrna_md_defaults_circ(flag)
    ccall((:vrna_md_defaults_circ, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_circ_get()
    ccall((:vrna_md_defaults_circ_get, libRNA), Cint, ())
end

function vrna_md_defaults_gquad(flag)
    ccall((:vrna_md_defaults_gquad, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_gquad_get()
    ccall((:vrna_md_defaults_gquad_get, libRNA), Cint, ())
end

function vrna_md_defaults_uniq_ML(flag)
    ccall((:vrna_md_defaults_uniq_ML, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_uniq_ML_get()
    ccall((:vrna_md_defaults_uniq_ML_get, libRNA), Cint, ())
end

function vrna_md_defaults_energy_set(e)
    ccall((:vrna_md_defaults_energy_set, libRNA), Cvoid, (Cint,), e)
end

function vrna_md_defaults_energy_set_get()
    ccall((:vrna_md_defaults_energy_set_get, libRNA), Cint, ())
end

function vrna_md_defaults_backtrack(flag)
    ccall((:vrna_md_defaults_backtrack, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_backtrack_get()
    ccall((:vrna_md_defaults_backtrack_get, libRNA), Cint, ())
end

function vrna_md_defaults_backtrack_type(t)
    ccall((:vrna_md_defaults_backtrack_type, libRNA), Cvoid, (Cchar,), t)
end

function vrna_md_defaults_backtrack_type_get()
    ccall((:vrna_md_defaults_backtrack_type_get, libRNA), Cchar, ())
end

function vrna_md_defaults_compute_bpp(flag)
    ccall((:vrna_md_defaults_compute_bpp, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_compute_bpp_get()
    ccall((:vrna_md_defaults_compute_bpp_get, libRNA), Cint, ())
end

function vrna_md_defaults_max_bp_span(span)
    ccall((:vrna_md_defaults_max_bp_span, libRNA), Cvoid, (Cint,), span)
end

function vrna_md_defaults_max_bp_span_get()
    ccall((:vrna_md_defaults_max_bp_span_get, libRNA), Cint, ())
end

function vrna_md_defaults_min_loop_size(size)
    ccall((:vrna_md_defaults_min_loop_size, libRNA), Cvoid, (Cint,), size)
end

function vrna_md_defaults_min_loop_size_get()
    ccall((:vrna_md_defaults_min_loop_size_get, libRNA), Cint, ())
end

function vrna_md_defaults_window_size(size)
    ccall((:vrna_md_defaults_window_size, libRNA), Cvoid, (Cint,), size)
end

function vrna_md_defaults_window_size_get()
    ccall((:vrna_md_defaults_window_size_get, libRNA), Cint, ())
end

function vrna_md_defaults_oldAliEn(flag)
    ccall((:vrna_md_defaults_oldAliEn, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_oldAliEn_get()
    ccall((:vrna_md_defaults_oldAliEn_get, libRNA), Cint, ())
end

function vrna_md_defaults_ribo(flag)
    ccall((:vrna_md_defaults_ribo, libRNA), Cvoid, (Cint,), flag)
end

function vrna_md_defaults_ribo_get()
    ccall((:vrna_md_defaults_ribo_get, libRNA), Cint, ())
end

function vrna_md_defaults_cv_fact(factor)
    ccall((:vrna_md_defaults_cv_fact, libRNA), Cvoid, (Cdouble,), factor)
end

function vrna_md_defaults_cv_fact_get()
    ccall((:vrna_md_defaults_cv_fact_get, libRNA), Cdouble, ())
end

function vrna_md_defaults_nc_fact(factor)
    ccall((:vrna_md_defaults_nc_fact, libRNA), Cvoid, (Cdouble,), factor)
end

function vrna_md_defaults_nc_fact_get()
    ccall((:vrna_md_defaults_nc_fact_get, libRNA), Cdouble, ())
end

function vrna_md_defaults_sfact(factor)
    ccall((:vrna_md_defaults_sfact, libRNA), Cvoid, (Cdouble,), factor)
end

function vrna_md_defaults_sfact_get()
    ccall((:vrna_md_defaults_sfact_get, libRNA), Cdouble, ())
end

function set_model_details(md)
    ccall((:set_model_details, libRNA), Cvoid, (Ptr{vrna_md_t},), md)
end

function option_string()
    ccall((:option_string, libRNA), Ptr{Cchar}, ())
end

struct vrna_param_s
    id::Cint
    stack::NTuple{8, NTuple{8, Cint}}
    hairpin::NTuple{31, Cint}
    bulge::NTuple{31, Cint}
    internal_loop::NTuple{31, Cint}
    mismatchExt::NTuple{8, NTuple{5, NTuple{5, Cint}}}
    mismatchI::NTuple{8, NTuple{5, NTuple{5, Cint}}}
    mismatch1nI::NTuple{8, NTuple{5, NTuple{5, Cint}}}
    mismatch23I::NTuple{8, NTuple{5, NTuple{5, Cint}}}
    mismatchH::NTuple{8, NTuple{5, NTuple{5, Cint}}}
    mismatchM::NTuple{8, NTuple{5, NTuple{5, Cint}}}
    dangle5::NTuple{8, NTuple{5, Cint}}
    dangle3::NTuple{8, NTuple{5, Cint}}
    int11::NTuple{8, NTuple{8, NTuple{5, NTuple{5, Cint}}}}
    int21::NTuple{8, NTuple{8, NTuple{5, NTuple{5, NTuple{5, Cint}}}}}
    int22::NTuple{8, NTuple{8, NTuple{5, NTuple{5, NTuple{5, NTuple{5, Cint}}}}}}
    ninio::NTuple{5, Cint}
    lxc::Cdouble
    MLbase::Cint
    MLintern::NTuple{8, Cint}
    MLclosing::Cint
    TerminalAU::Cint
    DuplexInit::Cint
    Tetraloop_E::NTuple{200, Cint}
    Tetraloops::NTuple{1401, Cchar}
    Triloop_E::NTuple{40, Cint}
    Triloops::NTuple{241, Cchar}
    Hexaloop_E::NTuple{40, Cint}
    Hexaloops::NTuple{1801, Cchar}
    TripleC::Cint
    MultipleCA::Cint
    MultipleCB::Cint
    gquad::NTuple{8, NTuple{46, Cint}}
    gquadLayerMismatch::Cint
    gquadLayerMismatchMax::Cint
    temperature::Cdouble
    model_details::vrna_md_t
    param_file::NTuple{256, Cchar}
end

const vrna_param_t = vrna_param_s

struct vrna_exp_param_s
    id::Cint
    expstack::NTuple{8, NTuple{8, Cdouble}}
    exphairpin::NTuple{31, Cdouble}
    expbulge::NTuple{31, Cdouble}
    expinternal::NTuple{31, Cdouble}
    expmismatchExt::NTuple{8, NTuple{5, NTuple{5, Cdouble}}}
    expmismatchI::NTuple{8, NTuple{5, NTuple{5, Cdouble}}}
    expmismatch23I::NTuple{8, NTuple{5, NTuple{5, Cdouble}}}
    expmismatch1nI::NTuple{8, NTuple{5, NTuple{5, Cdouble}}}
    expmismatchH::NTuple{8, NTuple{5, NTuple{5, Cdouble}}}
    expmismatchM::NTuple{8, NTuple{5, NTuple{5, Cdouble}}}
    expdangle5::NTuple{8, NTuple{5, Cdouble}}
    expdangle3::NTuple{8, NTuple{5, Cdouble}}
    expint11::NTuple{8, NTuple{8, NTuple{5, NTuple{5, Cdouble}}}}
    expint21::NTuple{8, NTuple{8, NTuple{5, NTuple{5, NTuple{5, Cdouble}}}}}
    expint22::NTuple{8, NTuple{8, NTuple{5, NTuple{5, NTuple{5, NTuple{5, Cdouble}}}}}}
    expninio::NTuple{5, NTuple{31, Cdouble}}
    lxc::Cdouble
    expMLbase::Cdouble
    expMLintern::NTuple{8, Cdouble}
    expMLclosing::Cdouble
    expTermAU::Cdouble
    expDuplexInit::Cdouble
    exptetra::NTuple{40, Cdouble}
    exptri::NTuple{40, Cdouble}
    exphex::NTuple{40, Cdouble}
    Tetraloops::NTuple{1401, Cchar}
    expTriloop::NTuple{40, Cdouble}
    Triloops::NTuple{241, Cchar}
    Hexaloops::NTuple{1801, Cchar}
    expTripleC::Cdouble
    expMultipleCA::Cdouble
    expMultipleCB::Cdouble
    expgquad::NTuple{8, NTuple{46, Cdouble}}
    expgquadLayerMismatch::Cdouble
    gquadLayerMismatchMax::Cint
    kT::Cdouble
    pf_scale::Cdouble
    temperature::Cdouble
    alpha::Cdouble
    model_details::vrna_md_t
    param_file::NTuple{256, Cchar}
end

const vrna_exp_param_t = vrna_exp_param_s

function vrna_params(md)
    ccall((:vrna_params, libRNA), Ptr{vrna_param_t}, (Ptr{vrna_md_t},), md)
end

function vrna_params_copy(par)
    ccall((:vrna_params_copy, libRNA), Ptr{vrna_param_t}, (Ptr{vrna_param_t},), par)
end

function vrna_exp_params(md)
    ccall((:vrna_exp_params, libRNA), Ptr{vrna_exp_param_t}, (Ptr{vrna_md_t},), md)
end

function vrna_exp_params_comparative(n_seq, md)
    ccall((:vrna_exp_params_comparative, libRNA), Ptr{vrna_exp_param_t}, (Cuint, Ptr{vrna_md_t}), n_seq, md)
end

function vrna_exp_params_copy(par)
    ccall((:vrna_exp_params_copy, libRNA), Ptr{vrna_exp_param_t}, (Ptr{vrna_exp_param_t},), par)
end

@cenum vrna_fc_type_e::UInt32 begin
    VRNA_FC_TYPE_SINGLE = 0
    VRNA_FC_TYPE_COMPARATIVE = 1
end

@cenum vrna_seq_type_e::UInt32 begin
    VRNA_SEQ_UNKNOWN = 0
    VRNA_SEQ_RNA = 1
    VRNA_SEQ_DNA = 2
end

struct vrna_sequence_s
    type::vrna_seq_type_e
    name::Ptr{Cchar}
    string::Ptr{Cchar}
    encoding::Ptr{Cshort}
    encoding5::Ptr{Cshort}
    encoding3::Ptr{Cshort}
    length::Cuint
end

const vrna_seq_t = vrna_sequence_s

struct vrna_alignment_s
    n_seq::Cuint
    sequences::Ptr{vrna_seq_t}
    gapfree_seq::Ptr{Ptr{Cchar}}
    gapfree_size::Ptr{Cuint}
    genome_size::Ptr{Culonglong}
    start::Ptr{Culonglong}
    orientation::Ptr{Cuchar}
    a2s::Ptr{Ptr{Cuint}}
end

const vrna_msa_t = vrna_alignment_s

@cenum vrna_hc_type_e::UInt32 begin
    VRNA_HC_DEFAULT = 0
    VRNA_HC_WINDOW = 1
end

mutable struct vrna_hc_depot_s end

const vrna_hc_depot_t = vrna_hc_depot_s

struct vrna_hc_s
    type::vrna_hc_type_e
    n::Cuint
    state::Cuchar
    mx::Ptr{Cuchar}
    matrix_local::Ptr{Ptr{Cuchar}}
    up_ext::Ptr{Cint}
    up_hp::Ptr{Cint}
    up_int::Ptr{Cint}
    up_ml::Ptr{Cint}
    f::Ptr{Cvoid}
    data::Ptr{Cvoid}
    free_data::Ptr{Cvoid}
    depot::Ptr{vrna_hc_depot_t}
end

const vrna_hc_t = vrna_hc_s

@cenum vrna_mx_type_e::UInt32 begin
    VRNA_MX_DEFAULT = 0
    VRNA_MX_WINDOW = 1
    VRNA_MX_2DFOLD = 2
end

struct vrna_mx_mfe_s
    type::vrna_mx_type_e
    length::Cuint
    strands::Cuint
    c::Ptr{Cint}
    f5::Ptr{Cint}
    f3::Ptr{Cint}
    fms5::Ptr{Ptr{Cint}}
    fms3::Ptr{Ptr{Cint}}
    fML::Ptr{Cint}
    fM1::Ptr{Cint}
    fM2::Ptr{Cint}
    ggg::Ptr{Cint}
    Fc::Cint
    FcH::Cint
    FcI::Cint
    FcM::Cint
    c_local::Ptr{Ptr{Cint}}
    f3_local::Ptr{Cint}
    fML_local::Ptr{Ptr{Cint}}
    ggg_local::Ptr{Ptr{Cint}}
    E_F5::Ptr{Ptr{Ptr{Cint}}}
    l_min_F5::Ptr{Ptr{Cint}}
    l_max_F5::Ptr{Ptr{Cint}}
    k_min_F5::Ptr{Cint}
    k_max_F5::Ptr{Cint}
    E_F3::Ptr{Ptr{Ptr{Cint}}}
    l_min_F3::Ptr{Ptr{Cint}}
    l_max_F3::Ptr{Ptr{Cint}}
    k_min_F3::Ptr{Cint}
    k_max_F3::Ptr{Cint}
    E_C::Ptr{Ptr{Ptr{Cint}}}
    l_min_C::Ptr{Ptr{Cint}}
    l_max_C::Ptr{Ptr{Cint}}
    k_min_C::Ptr{Cint}
    k_max_C::Ptr{Cint}
    E_M::Ptr{Ptr{Ptr{Cint}}}
    l_min_M::Ptr{Ptr{Cint}}
    l_max_M::Ptr{Ptr{Cint}}
    k_min_M::Ptr{Cint}
    k_max_M::Ptr{Cint}
    E_M1::Ptr{Ptr{Ptr{Cint}}}
    l_min_M1::Ptr{Ptr{Cint}}
    l_max_M1::Ptr{Ptr{Cint}}
    k_min_M1::Ptr{Cint}
    k_max_M1::Ptr{Cint}
    E_M2::Ptr{Ptr{Ptr{Cint}}}
    l_min_M2::Ptr{Ptr{Cint}}
    l_max_M2::Ptr{Ptr{Cint}}
    k_min_M2::Ptr{Cint}
    k_max_M2::Ptr{Cint}
    E_Fc::Ptr{Ptr{Cint}}
    l_min_Fc::Ptr{Cint}
    l_max_Fc::Ptr{Cint}
    k_min_Fc::Cint
    k_max_Fc::Cint
    E_FcH::Ptr{Ptr{Cint}}
    l_min_FcH::Ptr{Cint}
    l_max_FcH::Ptr{Cint}
    k_min_FcH::Cint
    k_max_FcH::Cint
    E_FcI::Ptr{Ptr{Cint}}
    l_min_FcI::Ptr{Cint}
    l_max_FcI::Ptr{Cint}
    k_min_FcI::Cint
    k_max_FcI::Cint
    E_FcM::Ptr{Ptr{Cint}}
    l_min_FcM::Ptr{Cint}
    l_max_FcM::Ptr{Cint}
    k_min_FcM::Cint
    k_max_FcM::Cint
    E_F5_rem::Ptr{Cint}
    E_F3_rem::Ptr{Cint}
    E_C_rem::Ptr{Cint}
    E_M_rem::Ptr{Cint}
    E_M1_rem::Ptr{Cint}
    E_M2_rem::Ptr{Cint}
    E_Fc_rem::Cint
    E_FcH_rem::Cint
    E_FcI_rem::Cint
    E_FcM_rem::Cint
end

const vrna_mx_mfe_t = vrna_mx_mfe_s

struct vrna_mx_pf_s
    type::vrna_mx_type_e
    length::Cuint
    scale::Ptr{FLT_OR_DBL}
    expMLbase::Ptr{FLT_OR_DBL}
    q::Ptr{FLT_OR_DBL}
    qb::Ptr{FLT_OR_DBL}
    qm::Ptr{FLT_OR_DBL}
    qm1::Ptr{FLT_OR_DBL}
    probs::Ptr{FLT_OR_DBL}
    q1k::Ptr{FLT_OR_DBL}
    qln::Ptr{FLT_OR_DBL}
    G::Ptr{FLT_OR_DBL}
    qo::FLT_OR_DBL
    qm2::Ptr{FLT_OR_DBL}
    qho::FLT_OR_DBL
    qio::FLT_OR_DBL
    qmo::FLT_OR_DBL
    q_local::Ptr{Ptr{FLT_OR_DBL}}
    qb_local::Ptr{Ptr{FLT_OR_DBL}}
    qm_local::Ptr{Ptr{FLT_OR_DBL}}
    pR::Ptr{Ptr{FLT_OR_DBL}}
    qm2_local::Ptr{Ptr{FLT_OR_DBL}}
    QI5::Ptr{Ptr{FLT_OR_DBL}}
    q2l::Ptr{Ptr{FLT_OR_DBL}}
    qmb::Ptr{Ptr{FLT_OR_DBL}}
    G_local::Ptr{Ptr{FLT_OR_DBL}}
    Q::Ptr{Ptr{Ptr{FLT_OR_DBL}}}
    l_min_Q::Ptr{Ptr{Cint}}
    l_max_Q::Ptr{Ptr{Cint}}
    k_min_Q::Ptr{Cint}
    k_max_Q::Ptr{Cint}
    Q_B::Ptr{Ptr{Ptr{FLT_OR_DBL}}}
    l_min_Q_B::Ptr{Ptr{Cint}}
    l_max_Q_B::Ptr{Ptr{Cint}}
    k_min_Q_B::Ptr{Cint}
    k_max_Q_B::Ptr{Cint}
    Q_M::Ptr{Ptr{Ptr{FLT_OR_DBL}}}
    l_min_Q_M::Ptr{Ptr{Cint}}
    l_max_Q_M::Ptr{Ptr{Cint}}
    k_min_Q_M::Ptr{Cint}
    k_max_Q_M::Ptr{Cint}
    Q_M1::Ptr{Ptr{Ptr{FLT_OR_DBL}}}
    l_min_Q_M1::Ptr{Ptr{Cint}}
    l_max_Q_M1::Ptr{Ptr{Cint}}
    k_min_Q_M1::Ptr{Cint}
    k_max_Q_M1::Ptr{Cint}
    Q_M2::Ptr{Ptr{Ptr{FLT_OR_DBL}}}
    l_min_Q_M2::Ptr{Ptr{Cint}}
    l_max_Q_M2::Ptr{Ptr{Cint}}
    k_min_Q_M2::Ptr{Cint}
    k_max_Q_M2::Ptr{Cint}
    Q_c::Ptr{Ptr{FLT_OR_DBL}}
    l_min_Q_c::Ptr{Cint}
    l_max_Q_c::Ptr{Cint}
    k_min_Q_c::Cint
    k_max_Q_c::Cint
    Q_cH::Ptr{Ptr{FLT_OR_DBL}}
    l_min_Q_cH::Ptr{Cint}
    l_max_Q_cH::Ptr{Cint}
    k_min_Q_cH::Cint
    k_max_Q_cH::Cint
    Q_cI::Ptr{Ptr{FLT_OR_DBL}}
    l_min_Q_cI::Ptr{Cint}
    l_max_Q_cI::Ptr{Cint}
    k_min_Q_cI::Cint
    k_max_Q_cI::Cint
    Q_cM::Ptr{Ptr{FLT_OR_DBL}}
    l_min_Q_cM::Ptr{Cint}
    l_max_Q_cM::Ptr{Cint}
    k_min_Q_cM::Cint
    k_max_Q_cM::Cint
    Q_rem::Ptr{FLT_OR_DBL}
    Q_B_rem::Ptr{FLT_OR_DBL}
    Q_M_rem::Ptr{FLT_OR_DBL}
    Q_M1_rem::Ptr{FLT_OR_DBL}
    Q_M2_rem::Ptr{FLT_OR_DBL}
    Q_c_rem::FLT_OR_DBL
    Q_cH_rem::FLT_OR_DBL
    Q_cI_rem::FLT_OR_DBL
    Q_cM_rem::FLT_OR_DBL
end

const vrna_mx_pf_t = vrna_mx_pf_s

struct vrna_structured_domains_s
    __placeholder::Cchar
end

const vrna_sd_t = vrna_structured_domains_s

struct vrna_unstructured_domain_s
    uniq_motif_count::Cint
    uniq_motif_size::Ptr{Cuint}
    motif_count::Cint
    motif::Ptr{Ptr{Cchar}}
    motif_name::Ptr{Ptr{Cchar}}
    motif_size::Ptr{Cuint}
    motif_en::Ptr{Cdouble}
    motif_type::Ptr{Cuint}
    prod_cb::Ptr{Cvoid}
    exp_prod_cb::Ptr{Cvoid}
    energy_cb::Ptr{Cvoid}
    exp_energy_cb::Ptr{Cvoid}
    data::Ptr{Cvoid}
    free_data::Ptr{Cvoid}
    probs_add::Ptr{Cvoid}
    probs_get::Ptr{Cvoid}
end

const vrna_ud_t = vrna_unstructured_domain_s

struct vrna_gr_aux_s
    cb_proc::Ptr{Cvoid}
    cb_aux_f::Ptr{Cvoid}
    cb_aux_c::Ptr{Cvoid}
    cb_aux_m::Ptr{Cvoid}
    cb_aux_m1::Ptr{Cvoid}
    cb_aux::Ptr{Cvoid}
    cb_aux_exp_f::Ptr{Cvoid}
    cb_aux_exp_c::Ptr{Cvoid}
    cb_aux_exp_m::Ptr{Cvoid}
    cb_aux_exp_m1::Ptr{Cvoid}
    cb_aux_exp::Ptr{Cvoid}
    data::Ptr{Cvoid}
    free_data::Ptr{Cvoid}
end

const vrna_gr_aux_t = vrna_gr_aux_s

@cenum vrna_sc_type_e::UInt32 begin
    VRNA_SC_DEFAULT = 0
    VRNA_SC_WINDOW = 1
end

struct vrna_sc_bp_storage_t
    interval_start::Cuint
    interval_end::Cuint
    e::Cint
end

struct vrna_sc_s
    type::vrna_sc_type_e
    n::Cuint
    state::Cuchar
    energy_up::Ptr{Ptr{Cint}}
    exp_energy_up::Ptr{Ptr{FLT_OR_DBL}}
    up_storage::Ptr{Cint}
    bp_storage::Ptr{Ptr{vrna_sc_bp_storage_t}}
    energy_bp::Ptr{Cint}
    exp_energy_bp::Ptr{FLT_OR_DBL}
    energy_bp_local::Ptr{Ptr{Cint}}
    exp_energy_bp_local::Ptr{Ptr{FLT_OR_DBL}}
    energy_stack::Ptr{Cint}
    exp_energy_stack::Ptr{FLT_OR_DBL}
    f::Ptr{Cvoid}
    bt::Ptr{Cvoid}
    exp_f::Ptr{Cvoid}
    data::Ptr{Cvoid}
    free_data::Ptr{Cvoid}
end

const vrna_sc_t = vrna_sc_s

struct vrna_fc_s
    type::vrna_fc_type_e
    length::Cuint
    cutpoint::Cint
    strand_number::Ptr{Cuint}
    strand_order::Ptr{Cuint}
    strand_order_uniq::Ptr{Cuint}
    strand_start::Ptr{Cuint}
    strand_end::Ptr{Cuint}
    strands::Cuint
    nucleotides::Ptr{vrna_seq_t}
    alignment::Ptr{vrna_msa_t}
    hc::Ptr{vrna_hc_t}
    matrices::Ptr{vrna_mx_mfe_t}
    exp_matrices::Ptr{vrna_mx_pf_t}
    params::Ptr{vrna_param_t}
    exp_params::Ptr{vrna_exp_param_t}
    iindx::Ptr{Cint}
    jindx::Ptr{Cint}
    stat_cb::Ptr{Cvoid}
    auxdata::Ptr{Cvoid}
    free_auxdata::Ptr{Cvoid}
    domains_struc::Ptr{vrna_sd_t}
    domains_up::Ptr{vrna_ud_t}
    aux_grammar::Ptr{vrna_gr_aux_t}
    sequence::Ptr{Cchar}
    sequence_encoding::Ptr{Cshort}
    sequence_encoding2::Ptr{Cshort}
    ptype::Ptr{Cchar}
    ptype_pf_compat::Ptr{Cchar}
    sc::Ptr{vrna_sc_t}
    sequences::Ptr{Ptr{Cchar}}
    n_seq::Cuint
    cons_seq::Ptr{Cchar}
    S_cons::Ptr{Cshort}
    S::Ptr{Ptr{Cshort}}
    S5::Ptr{Ptr{Cshort}}
    S3::Ptr{Ptr{Cshort}}
    Ss::Ptr{Ptr{Cchar}}
    a2s::Ptr{Ptr{Cuint}}
    pscore::Ptr{Cint}
    pscore_local::Ptr{Ptr{Cint}}
    pscore_pf_compat::Ptr{Cshort}
    scs::Ptr{Ptr{vrna_sc_t}}
    oldAliEn::Cint
    maxD1::Cuint
    maxD2::Cuint
    reference_pt1::Ptr{Cshort}
    reference_pt2::Ptr{Cshort}
    referenceBPs1::Ptr{Cuint}
    referenceBPs2::Ptr{Cuint}
    bpdist::Ptr{Cuint}
    mm1::Ptr{Cuint}
    mm2::Ptr{Cuint}
    window_size::Cint
    ptype_local::Ptr{Ptr{Cchar}}
end

const vrna_fold_compound_t = vrna_fc_s

function vrna_params_subst(vc, par)
    ccall((:vrna_params_subst, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{vrna_param_t}), vc, par)
end

function vrna_exp_params_subst(vc, params)
    ccall((:vrna_exp_params_subst, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{vrna_exp_param_t}), vc, params)
end

function vrna_exp_params_rescale(vc, mfe)
    ccall((:vrna_exp_params_rescale, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cdouble}), vc, mfe)
end

function vrna_params_reset(vc, md_p)
    ccall((:vrna_params_reset, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{vrna_md_t}), vc, md_p)
end

function vrna_exp_params_reset(vc, md_p)
    ccall((:vrna_exp_params_reset, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{vrna_md_t}), vc, md_p)
end

function vrna_params_prepare(vc, options)
    ccall((:vrna_params_prepare, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Cuint), vc, options)
end

const paramT = vrna_param_s

const pf_paramT = vrna_exp_param_s

function get_parameter_copy(par)
    ccall((:get_parameter_copy, libRNA), Ptr{vrna_param_t}, (Ptr{vrna_param_t},), par)
end

function get_scaled_pf_parameters()
    ccall((:get_scaled_pf_parameters, libRNA), Ptr{vrna_exp_param_t}, ())
end

function get_boltzmann_factors(temperature, betaScale, md, pf_scale)
    ccall((:get_boltzmann_factors, libRNA), Ptr{vrna_exp_param_t}, (Cdouble, Cdouble, vrna_md_t, Cdouble), temperature, betaScale, md, pf_scale)
end

function get_boltzmann_factor_copy(parameters)
    ccall((:get_boltzmann_factor_copy, libRNA), Ptr{vrna_exp_param_t}, (Ptr{vrna_exp_param_t},), parameters)
end

function get_scaled_alipf_parameters(n_seq)
    ccall((:get_scaled_alipf_parameters, libRNA), Ptr{vrna_exp_param_t}, (Cuint,), n_seq)
end

function get_boltzmann_factors_ali(n_seq, temperature, betaScale, md, pf_scale)
    ccall((:get_boltzmann_factors_ali, libRNA), Ptr{vrna_exp_param_t}, (Cuint, Cdouble, Cdouble, vrna_md_t, Cdouble), n_seq, temperature, betaScale, md, pf_scale)
end

function scale_parameters()
    ccall((:scale_parameters, libRNA), Ptr{vrna_param_t}, ())
end

function get_scaled_parameters(temperature, md)
    ccall((:get_scaled_parameters, libRNA), Ptr{vrna_param_t}, (Cdouble, vrna_md_t), temperature, md)
end

function copy_parameters()
    ccall((:copy_parameters, libRNA), Ptr{vrna_param_t}, ())
end

function set_parameters(dest)
    ccall((:set_parameters, libRNA), Ptr{vrna_param_t}, (Ptr{vrna_param_t},), dest)
end

function scale_pf_parameters()
    ccall((:scale_pf_parameters, libRNA), Ptr{vrna_exp_param_t}, ())
end

function copy_pf_param()
    ccall((:copy_pf_param, libRNA), Ptr{vrna_exp_param_t}, ())
end

function set_pf_param(dest)
    ccall((:set_pf_param, libRNA), Ptr{vrna_exp_param_t}, (Ptr{vrna_param_t},), dest)
end

function vrna_sequence(string, options)
    ccall((:vrna_sequence, libRNA), Ptr{vrna_seq_t}, (Ptr{Cchar}, Cuint), string, options)
end

function vrna_sequence_add(fc, string, options)
    ccall((:vrna_sequence_add, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Cuint), fc, string, options)
end

function vrna_sequence_remove(fc, i)
    ccall((:vrna_sequence_remove, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cuint), fc, i)
end

function vrna_sequence_remove_all(fc)
    ccall((:vrna_sequence_remove_all, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), fc)
end

function vrna_sequence_prepare(fc)
    ccall((:vrna_sequence_prepare, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), fc)
end

function vrna_sequence_order_update(fc, order)
    ccall((:vrna_sequence_order_update, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cuint}), fc, order)
end

function vrna_msa_add(fc, alignment, names, orientation, start, genome_size, options)
    ccall((:vrna_msa_add, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}, Ptr{Cuchar}, Ptr{Culonglong}, Ptr{Culonglong}, Cuint), fc, alignment, names, orientation, start, genome_size, options)
end

function vrna_mx_add(vc, type, options)
    ccall((:vrna_mx_add, libRNA), Cint, (Ptr{vrna_fold_compound_t}, vrna_mx_type_e, Cuint), vc, type, options)
end

function vrna_mx_mfe_add(vc, mx_type, options)
    ccall((:vrna_mx_mfe_add, libRNA), Cint, (Ptr{vrna_fold_compound_t}, vrna_mx_type_e, Cuint), vc, mx_type, options)
end

function vrna_mx_pf_add(vc, mx_type, options)
    ccall((:vrna_mx_pf_add, libRNA), Cint, (Ptr{vrna_fold_compound_t}, vrna_mx_type_e, Cuint), vc, mx_type, options)
end

function vrna_mx_prepare(vc, options)
    ccall((:vrna_mx_prepare, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cuint), vc, options)
end

function vrna_mx_mfe_free(vc)
    ccall((:vrna_mx_mfe_free, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), vc)
end

function vrna_mx_pf_free(vc)
    ccall((:vrna_mx_pf_free, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), vc)
end

struct vrna_hc_up_s
    position::Cint
    strand::Cint
    options::Cuchar
end

const vrna_hc_up_t = vrna_hc_up_s

function vrna_constraints_add(vc, constraint, options)
    ccall((:vrna_constraints_add, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Cuint), vc, constraint, options)
end

# typedef unsigned char ( vrna_callback_hc_evaluate ) ( int i , int j , int k , int l , unsigned char d , void * data )
const vrna_callback_hc_evaluate = Cvoid

function vrna_message_constraint_options(option)
    ccall((:vrna_message_constraint_options, libRNA), Cvoid, (Cuint,), option)
end

function vrna_message_constraint_options_all()
    ccall((:vrna_message_constraint_options_all, libRNA), Cvoid, ())
end

function vrna_hc_init(vc)
    ccall((:vrna_hc_init, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), vc)
end

function vrna_hc_init_window(vc)
    ccall((:vrna_hc_init_window, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), vc)
end

function vrna_hc_prepare(fc, options)
    ccall((:vrna_hc_prepare, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cuint), fc, options)
end

function vrna_hc_update(fc, i, options)
    ccall((:vrna_hc_update, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Cuint, Cuint), fc, i, options)
end

function vrna_hc_add_up(vc, i, option)
    ccall((:vrna_hc_add_up, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Cint, Cuchar), vc, i, option)
end

function vrna_hc_add_up_strand(fc, i, strand, option)
    ccall((:vrna_hc_add_up_strand, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Cuchar), fc, i, strand, option)
end

function vrna_hc_add_up_batch(vc, constraints)
    ccall((:vrna_hc_add_up_batch, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{vrna_hc_up_t}), vc, constraints)
end

function vrna_hc_add_up_strand_batch(fc, constraints)
    ccall((:vrna_hc_add_up_strand_batch, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{vrna_hc_up_t}), fc, constraints)
end

function vrna_hc_add_bp(vc, i, j, option)
    ccall((:vrna_hc_add_bp, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cint, Cint, Cuchar), vc, i, j, option)
end

function vrna_hc_add_bp_strand(fc, i, strand_i, j, strand_j, option)
    ccall((:vrna_hc_add_bp_strand, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Cuint, Cuint, Cuchar), fc, i, strand_i, j, strand_j, option)
end

function vrna_hc_add_bp_nonspecific(vc, i, d, option)
    ccall((:vrna_hc_add_bp_nonspecific, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Cint, Cint, Cuchar), vc, i, d, option)
end

function vrna_hc_free(hc)
    ccall((:vrna_hc_free, libRNA), Cvoid, (Ptr{vrna_hc_t},), hc)
end

function vrna_hc_add_f(vc, f)
    ccall((:vrna_hc_add_f, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), vc, f)
end

function vrna_hc_add_data(vc, data, f)
    ccall((:vrna_hc_add_data, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}, Ptr{Cvoid}), vc, data, f)
end

function vrna_hc_add_from_db(vc, constraint, options)
    ccall((:vrna_hc_add_from_db, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Cuint), vc, constraint, options)
end

function print_tty_constraint(option)
    ccall((:print_tty_constraint, libRNA), Cvoid, (Cuint,), option)
end

function print_tty_constraint_full()
    ccall((:print_tty_constraint_full, libRNA), Cvoid, ())
end

function constrain_ptypes(constraint, length, ptype, BP, min_loop_size, idx_type)
    ccall((:constrain_ptypes, libRNA), Cvoid, (Ptr{Cchar}, Cuint, Ptr{Cchar}, Ptr{Cint}, Cint, Cuint), constraint, length, ptype, BP, min_loop_size, idx_type)
end

# typedef int ( vrna_callback_sc_energy ) ( int i , int j , int k , int l , unsigned char d , void * data )
const vrna_callback_sc_energy = Cvoid

# typedef FLT_OR_DBL ( vrna_callback_sc_exp_energy ) ( int i , int j , int k , int l , unsigned char d , void * data )
const vrna_callback_sc_exp_energy = Cvoid

# typedef vrna_basepair_t * ( vrna_callback_sc_backtrack ) ( int i , int j , int k , int l , unsigned char d , void * data )
const vrna_callback_sc_backtrack = Cvoid

function vrna_sc_init(vc)
    ccall((:vrna_sc_init, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), vc)
end

function vrna_sc_prepare(vc, options)
    ccall((:vrna_sc_prepare, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Cuint), vc, options)
end

function vrna_sc_update(vc, i, options)
    ccall((:vrna_sc_update, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cuint, Cuint), vc, i, options)
end

function vrna_sc_set_bp(vc, constraints, options)
    ccall((:vrna_sc_set_bp, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Ptr{FLT_OR_DBL}}, Cuint), vc, constraints, options)
end

function vrna_sc_add_bp(vc, i, j, energy, options)
    ccall((:vrna_sc_add_bp, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cint, Cint, FLT_OR_DBL, Cuint), vc, i, j, energy, options)
end

function vrna_sc_set_up(vc, constraints, options)
    ccall((:vrna_sc_set_up, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{FLT_OR_DBL}, Cuint), vc, constraints, options)
end

function vrna_sc_add_up(vc, i, energy, options)
    ccall((:vrna_sc_add_up, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cint, FLT_OR_DBL, Cuint), vc, i, energy, options)
end

function vrna_sc_set_stack(vc, constraints, options)
    ccall((:vrna_sc_set_stack, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{FLT_OR_DBL}, Cuint), vc, constraints, options)
end

function vrna_sc_set_stack_comparative(fc, constraints, options)
    ccall((:vrna_sc_set_stack_comparative, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Ptr{FLT_OR_DBL}}, Cuint), fc, constraints, options)
end

function vrna_sc_add_stack(vc, i, energy, options)
    ccall((:vrna_sc_add_stack, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cint, FLT_OR_DBL, Cuint), vc, i, energy, options)
end

function vrna_sc_add_stack_comparative(fc, i, energies, options)
    ccall((:vrna_sc_add_stack_comparative, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cint, Ptr{FLT_OR_DBL}, Cuint), fc, i, energies, options)
end

function vrna_sc_remove(vc)
    ccall((:vrna_sc_remove, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), vc)
end

function vrna_sc_free(sc)
    ccall((:vrna_sc_free, libRNA), Cvoid, (Ptr{vrna_sc_t},), sc)
end

function vrna_sc_add_data(vc, data, free_data)
    ccall((:vrna_sc_add_data, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}, Ptr{Cvoid}), vc, data, free_data)
end

function vrna_sc_add_data_comparative(vc, data, free_data)
    ccall((:vrna_sc_add_data_comparative, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}), vc, data, free_data)
end

function vrna_sc_add_f(vc, f)
    ccall((:vrna_sc_add_f, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), vc, f)
end

function vrna_sc_add_f_comparative(vc, f)
    ccall((:vrna_sc_add_f_comparative, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Ptr{Cvoid}}), vc, f)
end

function vrna_sc_add_bt(vc, f)
    ccall((:vrna_sc_add_bt, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), vc, f)
end

function vrna_sc_add_exp_f(vc, exp_f)
    ccall((:vrna_sc_add_exp_f, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), vc, exp_f)
end

function vrna_sc_add_exp_f_comparative(vc, exp_f)
    ccall((:vrna_sc_add_exp_f_comparative, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Ptr{Cvoid}}), vc, exp_f)
end

# typedef int ( vrna_callback_gr_rule ) ( vrna_fold_compound_t * vc , int i , int j , void * data )
const vrna_callback_gr_rule = Cvoid

# typedef void ( vrna_callback_gr_rule_aux ) ( vrna_fold_compound_t * vc , int i , int j , void * data )
const vrna_callback_gr_rule_aux = Cvoid

# typedef FLT_OR_DBL ( vrna_callback_gr_rule_exp ) ( vrna_fold_compound_t * vc , int i , int j , void * data )
const vrna_callback_gr_rule_exp = Cvoid

# typedef void ( vrna_callback_gr_rule_aux_exp ) ( vrna_fold_compound_t * vc , int i , int j , void * data )
const vrna_callback_gr_rule_aux_exp = Cvoid

# typedef void ( vrna_callback_gr_cond ) ( vrna_fold_compound_t * fc , unsigned char stage , void * data )
const vrna_callback_gr_cond = Cvoid

# typedef void ( vrna_callback_gr_free_data ) ( void * data )
const vrna_callback_gr_free_data = Cvoid

function vrna_gr_set_aux_f(fc, cb)
    ccall((:vrna_gr_set_aux_f, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux_exp_f(fc, cb)
    ccall((:vrna_gr_set_aux_exp_f, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux_c(fc, cb)
    ccall((:vrna_gr_set_aux_c, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux_exp_c(fc, cb)
    ccall((:vrna_gr_set_aux_exp_c, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux_m(fc, cb)
    ccall((:vrna_gr_set_aux_m, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux_exp_m(fc, cb)
    ccall((:vrna_gr_set_aux_exp_m, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux_m1(fc, cb)
    ccall((:vrna_gr_set_aux_m1, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux_exp_m1(fc, cb)
    ccall((:vrna_gr_set_aux_exp_m1, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux(fc, cb)
    ccall((:vrna_gr_set_aux, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_aux_exp(fc, cb)
    ccall((:vrna_gr_set_aux_exp, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_set_data(fc, data, free_data)
    ccall((:vrna_gr_set_data, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}, Ptr{Cvoid}), fc, data, free_data)
end

function vrna_gr_set_cond(fc, cb)
    ccall((:vrna_gr_set_cond, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, cb)
end

function vrna_gr_reset(fc)
    ccall((:vrna_gr_reset, libRNA), Cint, (Ptr{vrna_fold_compound_t},), fc)
end

struct vrna_unstructured_domain_motif_s
    start::Cint
    number::Cint
end

const vrna_ud_motif_t = vrna_unstructured_domain_motif_s

struct vrna_hx_s
    start::Cuint
    _end::Cuint
    length::Cuint
    up5::Cuint
    up3::Cuint
end

const vrna_hx_t = vrna_hx_s

const vrna_ep_t = vrna_elem_prob_s

function vrna_db_pack(struc)
    ccall((:vrna_db_pack, libRNA), Ptr{Cchar}, (Ptr{Cchar},), struc)
end

function vrna_db_unpack(packed)
    ccall((:vrna_db_unpack, libRNA), Ptr{Cchar}, (Ptr{Cchar},), packed)
end

function vrna_db_flatten(structure, options)
    ccall((:vrna_db_flatten, libRNA), Cvoid, (Ptr{Cchar}, Cuint), structure, options)
end

function vrna_db_flatten_to(string, target, options)
    ccall((:vrna_db_flatten_to, libRNA), Cvoid, (Ptr{Cchar}, Ptr{Cchar}, Cuint), string, target, options)
end

function vrna_db_from_ptable(pt)
    ccall((:vrna_db_from_ptable, libRNA), Ptr{Cchar}, (Ptr{Cshort},), pt)
end

function vrna_db_from_plist(pairs, n)
    ccall((:vrna_db_from_plist, libRNA), Ptr{Cchar}, (Ptr{vrna_ep_t}, Cuint), pairs, n)
end

function vrna_db_to_element_string(structure)
    ccall((:vrna_db_to_element_string, libRNA), Ptr{Cchar}, (Ptr{Cchar},), structure)
end

function vrna_db_pk_remove(structure, options)
    ccall((:vrna_db_pk_remove, libRNA), Ptr{Cchar}, (Ptr{Cchar}, Cuint), structure, options)
end

function vrna_ptable(structure)
    ccall((:vrna_ptable, libRNA), Ptr{Cshort}, (Ptr{Cchar},), structure)
end

function vrna_ptable_from_string(string, options)
    ccall((:vrna_ptable_from_string, libRNA), Ptr{Cshort}, (Ptr{Cchar}, Cuint), string, options)
end

function vrna_pt_pk_get(structure)
    ccall((:vrna_pt_pk_get, libRNA), Ptr{Cshort}, (Ptr{Cchar},), structure)
end

function vrna_ptable_copy(pt)
    ccall((:vrna_ptable_copy, libRNA), Ptr{Cshort}, (Ptr{Cshort},), pt)
end

function vrna_pt_ali_get(structure)
    ccall((:vrna_pt_ali_get, libRNA), Ptr{Cshort}, (Ptr{Cchar},), structure)
end

function vrna_pt_snoop_get(structure)
    ccall((:vrna_pt_snoop_get, libRNA), Ptr{Cshort}, (Ptr{Cchar},), structure)
end

function vrna_pt_pk_remove(ptable, options)
    ccall((:vrna_pt_pk_remove, libRNA), Ptr{Cshort}, (Ptr{Cshort}, Cuint), ptable, options)
end

function vrna_plist(struc, pr)
    ccall((:vrna_plist, libRNA), Ptr{vrna_ep_t}, (Ptr{Cchar}, Cfloat), struc, pr)
end

function vrna_plist_from_probs(vc, cut_off)
    ccall((:vrna_plist_from_probs, libRNA), Ptr{vrna_ep_t}, (Ptr{vrna_fold_compound_t}, Cdouble), vc, cut_off)
end

function vrna_db_from_WUSS(wuss)
    ccall((:vrna_db_from_WUSS, libRNA), Ptr{Cchar}, (Ptr{Cchar},), wuss)
end

function vrna_abstract_shapes(structure, level)
    ccall((:vrna_abstract_shapes, libRNA), Ptr{Cchar}, (Ptr{Cchar}, Cuint), structure, level)
end

function vrna_abstract_shapes_pt(pt, level)
    ccall((:vrna_abstract_shapes_pt, libRNA), Ptr{Cchar}, (Ptr{Cshort}, Cuint), pt, level)
end

function vrna_hx_from_ptable(pt)
    ccall((:vrna_hx_from_ptable, libRNA), Ptr{vrna_hx_t}, (Ptr{Cshort},), pt)
end

function vrna_hx_merge(list, maxdist)
    ccall((:vrna_hx_merge, libRNA), Ptr{vrna_hx_t}, (Ptr{vrna_hx_t}, Cint), list, maxdist)
end

function vrna_loopidx_from_ptable(pt)
    ccall((:vrna_loopidx_from_ptable, libRNA), Ptr{Cint}, (Ptr{Cshort},), pt)
end

function vrna_bp_distance_pt(pt1, pt2)
    ccall((:vrna_bp_distance_pt, libRNA), Cint, (Ptr{Cshort}, Ptr{Cshort}), pt1, pt2)
end

function vrna_bp_distance(str1, str2)
    ccall((:vrna_bp_distance, libRNA), Cint, (Ptr{Cchar}, Ptr{Cchar}), str1, str2)
end

function vrna_dist_mountain(str1, str2, p)
    ccall((:vrna_dist_mountain, libRNA), Cdouble, (Ptr{Cchar}, Ptr{Cchar}, Cuint), str1, str2, p)
end

function vrna_refBPcnt_matrix(reference_pt, turn)
    ccall((:vrna_refBPcnt_matrix, libRNA), Ptr{Cuint}, (Ptr{Cshort}, Cuint), reference_pt, turn)
end

function vrna_refBPdist_matrix(pt1, pt2, turn)
    ccall((:vrna_refBPdist_matrix, libRNA), Ptr{Cuint}, (Ptr{Cshort}, Ptr{Cshort}, Cuint), pt1, pt2, turn)
end

function vrna_db_from_probs(pr, length)
    ccall((:vrna_db_from_probs, libRNA), Ptr{Cchar}, (Ptr{FLT_OR_DBL}, Cuint), pr, length)
end

function vrna_bpp_symbol(x)
    ccall((:vrna_bpp_symbol, libRNA), Cchar, (Ptr{Cfloat},), x)
end

function vrna_db_from_bp_stack(bp, length)
    ccall((:vrna_db_from_bp_stack, libRNA), Ptr{Cchar}, (Ptr{vrna_bp_stack_t}, Cuint), bp, length)
end

function vrna_letter_structure(structure, bp, length)
    ccall((:vrna_letter_structure, libRNA), Cvoid, (Ptr{Cchar}, Ptr{vrna_bp_stack_t}, Cuint), structure, bp, length)
end

function vrna_db_to_tree_string(structure, type)
    ccall((:vrna_db_to_tree_string, libRNA), Ptr{Cchar}, (Ptr{Cchar}, Cuint), structure, type)
end

function vrna_tree_string_unweight(structure)
    ccall((:vrna_tree_string_unweight, libRNA), Ptr{Cchar}, (Ptr{Cchar},), structure)
end

function vrna_tree_string_to_db(tree)
    ccall((:vrna_tree_string_to_db, libRNA), Ptr{Cchar}, (Ptr{Cchar},), tree)
end

function assign_plist_from_db(pl, struc, pr)
    ccall((:assign_plist_from_db, libRNA), Cvoid, (Ptr{Ptr{vrna_ep_t}}, Ptr{Cchar}, Cfloat), pl, struc, pr)
end

function pack_structure(struc)
    ccall((:pack_structure, libRNA), Ptr{Cchar}, (Ptr{Cchar},), struc)
end

function unpack_structure(packed)
    ccall((:unpack_structure, libRNA), Ptr{Cchar}, (Ptr{Cchar},), packed)
end

function make_pair_table(structure)
    ccall((:make_pair_table, libRNA), Ptr{Cshort}, (Ptr{Cchar},), structure)
end

function make_pair_table_pk(structure)
    ccall((:make_pair_table_pk, libRNA), Ptr{Cshort}, (Ptr{Cchar},), structure)
end

function copy_pair_table(pt)
    ccall((:copy_pair_table, libRNA), Ptr{Cshort}, (Ptr{Cshort},), pt)
end

function alimake_pair_table(structure)
    ccall((:alimake_pair_table, libRNA), Ptr{Cshort}, (Ptr{Cchar},), structure)
end

function make_pair_table_snoop(structure)
    ccall((:make_pair_table_snoop, libRNA), Ptr{Cshort}, (Ptr{Cchar},), structure)
end

function make_loop_index_pt(pt)
    ccall((:make_loop_index_pt, libRNA), Ptr{Cint}, (Ptr{Cshort},), pt)
end

function bp_distance(str1, str2)
    ccall((:bp_distance, libRNA), Cint, (Ptr{Cchar}, Ptr{Cchar}), str1, str2)
end

function make_referenceBP_array(reference_pt, turn)
    ccall((:make_referenceBP_array, libRNA), Ptr{Cuint}, (Ptr{Cshort}, Cuint), reference_pt, turn)
end

function compute_BPdifferences(pt1, pt2, turn)
    ccall((:compute_BPdifferences, libRNA), Ptr{Cuint}, (Ptr{Cshort}, Ptr{Cshort}, Cuint), pt1, pt2, turn)
end

function assign_plist_from_pr(pl, probs, length, cutoff)
    ccall((:assign_plist_from_pr, libRNA), Cvoid, (Ptr{Ptr{vrna_ep_t}}, Ptr{FLT_OR_DBL}, Cint, Cdouble), pl, probs, length, cutoff)
end

function parenthesis_structure(structure, bp, length)
    ccall((:parenthesis_structure, libRNA), Cvoid, (Ptr{Cchar}, Ptr{vrna_bp_stack_t}, Cint), structure, bp, length)
end

function parenthesis_zuker(structure, bp, length)
    ccall((:parenthesis_zuker, libRNA), Cvoid, (Ptr{Cchar}, Ptr{vrna_bp_stack_t}, Cint), structure, bp, length)
end

function letter_structure(structure, bp, length)
    ccall((:letter_structure, libRNA), Cvoid, (Ptr{Cchar}, Ptr{vrna_bp_stack_t}, Cint), structure, bp, length)
end

function bppm_to_structure(structure, pr, length)
    ccall((:bppm_to_structure, libRNA), Cvoid, (Ptr{Cchar}, Ptr{FLT_OR_DBL}, Cuint), structure, pr, length)
end

function bppm_symbol(x)
    ccall((:bppm_symbol, libRNA), Cchar, (Ptr{Cfloat},), x)
end

# typedef int ( vrna_callback_ud_energy ) ( vrna_fold_compound_t * vc , int i , int j , unsigned int loop_type , void * data )
const vrna_callback_ud_energy = Cvoid

# typedef FLT_OR_DBL ( vrna_callback_ud_exp_energy ) ( vrna_fold_compound_t * vc , int i , int j , unsigned int loop_type , void * data )
const vrna_callback_ud_exp_energy = Cvoid

# typedef void ( vrna_callback_ud_production ) ( vrna_fold_compound_t * vc , void * data )
const vrna_callback_ud_production = Cvoid

# typedef void ( vrna_callback_ud_exp_production ) ( vrna_fold_compound_t * vc , void * data )
const vrna_callback_ud_exp_production = Cvoid

# typedef void ( vrna_callback_ud_probs_add ) ( vrna_fold_compound_t * vc , int i , int j , unsigned int loop_type , FLT_OR_DBL exp_energy , void * data )
const vrna_callback_ud_probs_add = Cvoid

# typedef FLT_OR_DBL ( vrna_callback_ud_probs_get ) ( vrna_fold_compound_t * vc , int i , int j , unsigned int loop_type , int motif , void * data )
const vrna_callback_ud_probs_get = Cvoid

function vrna_ud_motifs_centroid(fc, structure)
    ccall((:vrna_ud_motifs_centroid, libRNA), Ptr{vrna_ud_motif_t}, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), fc, structure)
end

function vrna_ud_motifs_MEA(fc, structure, probability_list)
    ccall((:vrna_ud_motifs_MEA, libRNA), Ptr{vrna_ud_motif_t}, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Ptr{vrna_ep_t}), fc, structure, probability_list)
end

function vrna_ud_motifs_MFE(fc, structure)
    ccall((:vrna_ud_motifs_MFE, libRNA), Ptr{vrna_ud_motif_t}, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), fc, structure)
end

function vrna_ud_add_motif(vc, motif, motif_en, motif_name, loop_type)
    ccall((:vrna_ud_add_motif, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Cdouble, Ptr{Cchar}, Cuint), vc, motif, motif_en, motif_name, loop_type)
end

function vrna_ud_get_motif_size_at(vc, i, loop_type)
    ccall((:vrna_ud_get_motif_size_at, libRNA), Ptr{Cint}, (Ptr{vrna_fold_compound_t}, Cint, Cuint), vc, i, loop_type)
end

function vrna_ud_get_motifs_at(vc, i, loop_type)
    ccall((:vrna_ud_get_motifs_at, libRNA), Ptr{Cint}, (Ptr{vrna_fold_compound_t}, Cint, Cuint), vc, i, loop_type)
end

function vrna_ud_detect_motifs(vc, structure)
    ccall((:vrna_ud_detect_motifs, libRNA), Ptr{vrna_ud_motif_t}, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), vc, structure)
end

function vrna_ud_remove(vc)
    ccall((:vrna_ud_remove, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), vc)
end

function vrna_ud_set_data(vc, data, free_cb)
    ccall((:vrna_ud_set_data, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}, Ptr{Cvoid}), vc, data, free_cb)
end

function vrna_ud_set_prod_rule_cb(vc, pre_cb, e_cb)
    ccall((:vrna_ud_set_prod_rule_cb, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}, Ptr{Cvoid}), vc, pre_cb, e_cb)
end

function vrna_ud_set_exp_prod_rule_cb(vc, pre_cb, exp_e_cb)
    ccall((:vrna_ud_set_exp_prod_rule_cb, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}, Ptr{Cvoid}), vc, pre_cb, exp_e_cb)
end

function vrna_ud_set_prob_cb(vc, setter, getter)
    ccall((:vrna_ud_set_prob_cb, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}, Ptr{Cvoid}), vc, setter, getter)
end

struct pu_contrib
    H::Ptr{Ptr{Cdouble}}
    I::Ptr{Ptr{Cdouble}}
    M::Ptr{Ptr{Cdouble}}
    E::Ptr{Ptr{Cdouble}}
    length::Cint
    w::Cint
end

struct interact
    Pi::Ptr{Cdouble}
    Gi::Ptr{Cdouble}
    Gikjl::Cdouble
    Gikjl_wo::Cdouble
    i::Cint
    k::Cint
    j::Cint
    l::Cint
    length::Cint
end

struct pu_out
    len::Cint
    u_vals::Cint
    contribs::Cint
    header::Ptr{Ptr{Cchar}}
    u_values::Ptr{Ptr{Cdouble}}
end

struct constrain
    indx::Ptr{Cint}
    ptype::Ptr{Cchar}
end

struct duplexT
    i::Cint
    j::Cint
    _end::Cint
    structure::Ptr{Cchar}
    energy::Cdouble
    energy_backtrack::Cdouble
    opening_backtrack_x::Cdouble
    opening_backtrack_y::Cdouble
    offset::Cint
    dG1::Cdouble
    dG2::Cdouble
    ddG::Cdouble
    tb::Cint
    te::Cint
    qb::Cint
    qe::Cint
end

struct node
    k::Cint
    energy::Cint
    next::Ptr{node}
end

const folden = node

struct snoopT
    i::Cint
    j::Cint
    u::Cint
    structure::Ptr{Cchar}
    energy::Cfloat
    Duplex_El::Cfloat
    Duplex_Er::Cfloat
    Loop_E::Cfloat
    Loop_D::Cfloat
    pscd::Cfloat
    psct::Cfloat
    pscg::Cfloat
    Duplex_Ol::Cfloat
    Duplex_Or::Cfloat
    Duplex_Ot::Cfloat
    fullStemEnergy::Cfloat
end

struct dupVar
    i::Cint
    j::Cint
    _end::Cint
    pk_helix::Ptr{Cchar}
    structure::Ptr{Cchar}
    energy::Cdouble
    offset::Cint
    dG1::Cdouble
    dG2::Cdouble
    ddG::Cdouble
    tb::Cint
    te::Cint
    qb::Cint
    qe::Cint
    inactive::Cint
    processed::Cint
end

function vrna_pairing_probs(vc, structure)
    ccall((:vrna_pairing_probs, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), vc, structure)
end

function vrna_mean_bp_distance_pr(length, pr)
    ccall((:vrna_mean_bp_distance_pr, libRNA), Cdouble, (Cint, Ptr{FLT_OR_DBL}), length, pr)
end

function vrna_mean_bp_distance(vc)
    ccall((:vrna_mean_bp_distance, libRNA), Cdouble, (Ptr{vrna_fold_compound_t},), vc)
end

function vrna_ensemble_defect_pt(fc, pt)
    ccall((:vrna_ensemble_defect_pt, libRNA), Cdouble, (Ptr{vrna_fold_compound_t}, Ptr{Cshort}), fc, pt)
end

function vrna_ensemble_defect(fc, structure)
    ccall((:vrna_ensemble_defect, libRNA), Cdouble, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), fc, structure)
end

function vrna_positional_entropy(fc)
    ccall((:vrna_positional_entropy, libRNA), Ptr{Cdouble}, (Ptr{vrna_fold_compound_t},), fc)
end

function vrna_stack_prob(vc, cutoff)
    ccall((:vrna_stack_prob, libRNA), Ptr{vrna_ep_t}, (Ptr{vrna_fold_compound_t}, Cdouble), vc, cutoff)
end

function vrna_pf_dimer_probs(FAB, FA, FB, prAB, prA, prB, Alength, exp_params)
    ccall((:vrna_pf_dimer_probs, libRNA), Cvoid, (Cdouble, Cdouble, Cdouble, Ptr{vrna_ep_t}, Ptr{vrna_ep_t}, Ptr{vrna_ep_t}, Cint, Ptr{vrna_exp_param_t}), FAB, FA, FB, prAB, prA, prB, Alength, exp_params)
end

function vrna_pr_structure(fc, structure)
    ccall((:vrna_pr_structure, libRNA), Cdouble, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), fc, structure)
end

function vrna_pr_energy(vc, e)
    ccall((:vrna_pr_energy, libRNA), Cdouble, (Ptr{vrna_fold_compound_t}, Cdouble), vc, e)
end

mutable struct vrna_cstr_s end

const vrna_cstr_t = Ptr{vrna_cstr_s}

function vrna_cstr(size, output)
    ccall((:vrna_cstr, libRNA), vrna_cstr_t, (Csize_t, Ptr{Libc.FILE}), size, output)
end

function vrna_cstr_discard(buf)
    ccall((:vrna_cstr_discard, libRNA), Cvoid, (Ptr{vrna_cstr_s},), buf)
end

function vrna_cstr_free(buf)
    ccall((:vrna_cstr_free, libRNA), Cvoid, (vrna_cstr_t,), buf)
end

function vrna_cstr_close(buf)
    ccall((:vrna_cstr_close, libRNA), Cvoid, (vrna_cstr_t,), buf)
end

function vrna_cstr_fflush(buf)
    ccall((:vrna_cstr_fflush, libRNA), Cvoid, (Ptr{vrna_cstr_s},), buf)
end

function vrna_cstr_string(buf)
    ccall((:vrna_cstr_string, libRNA), Ptr{Cchar}, (vrna_cstr_t,), buf)
end

function vrna_cstr_print_fasta_header(buf, head)
    ccall((:vrna_cstr_print_fasta_header, libRNA), Cvoid, (vrna_cstr_t, Ptr{Cchar}), buf, head)
end

function vrna_cstr_print_eval_sd_corr(buf)
    ccall((:vrna_cstr_print_eval_sd_corr, libRNA), Cvoid, (Ptr{vrna_cstr_s},), buf)
end

function vrna_cstr_print_eval_ext_loop(buf, energy)
    ccall((:vrna_cstr_print_eval_ext_loop, libRNA), Cvoid, (Ptr{vrna_cstr_s}, Cint), buf, energy)
end

function vrna_cstr_print_eval_hp_loop(buf, i, j, si, sj, energy)
    ccall((:vrna_cstr_print_eval_hp_loop, libRNA), Cvoid, (Ptr{vrna_cstr_s}, Cint, Cint, Cchar, Cchar, Cint), buf, i, j, si, sj, energy)
end

function vrna_cstr_print_eval_hp_loop_revert(buf, i, j, si, sj, energy)
    ccall((:vrna_cstr_print_eval_hp_loop_revert, libRNA), Cvoid, (Ptr{vrna_cstr_s}, Cint, Cint, Cchar, Cchar, Cint), buf, i, j, si, sj, energy)
end

function vrna_cstr_print_eval_int_loop(buf, i, j, si, sj, k, l, sk, sl, energy)
    ccall((:vrna_cstr_print_eval_int_loop, libRNA), Cvoid, (Ptr{vrna_cstr_s}, Cint, Cint, Cchar, Cchar, Cint, Cint, Cchar, Cchar, Cint), buf, i, j, si, sj, k, l, sk, sl, energy)
end

function vrna_cstr_print_eval_int_loop_revert(buf, i, j, si, sj, k, l, sk, sl, energy)
    ccall((:vrna_cstr_print_eval_int_loop_revert, libRNA), Cvoid, (Ptr{vrna_cstr_s}, Cint, Cint, Cchar, Cchar, Cint, Cint, Cchar, Cchar, Cint), buf, i, j, si, sj, k, l, sk, sl, energy)
end

function vrna_cstr_print_eval_mb_loop(buf, i, j, si, sj, energy)
    ccall((:vrna_cstr_print_eval_mb_loop, libRNA), Cvoid, (Ptr{vrna_cstr_s}, Cint, Cint, Cchar, Cchar, Cint), buf, i, j, si, sj, energy)
end

function vrna_cstr_print_eval_mb_loop_revert(buf, i, j, si, sj, energy)
    ccall((:vrna_cstr_print_eval_mb_loop_revert, libRNA), Cvoid, (Ptr{vrna_cstr_s}, Cint, Cint, Cchar, Cchar, Cint), buf, i, j, si, sj, energy)
end

function vrna_cstr_print_eval_gquad(buf, i, L, l, energy)
    ccall((:vrna_cstr_print_eval_gquad, libRNA), Cvoid, (Ptr{vrna_cstr_s}, Cint, Cint, Ptr{Cint}, Cint), buf, i, L, l, energy)
end

struct vrna_move_s
    pos_5::Cint
    pos_3::Cint
    next::Ptr{Cvoid} # next::Ptr{vrna_move_t}
end

function Base.getproperty(x::vrna_move_s, f::Symbol)
    f === :next && return Ptr{vrna_move_t}(getfield(x, f))
    return getfield(x, f)
end

const vrna_move_t = vrna_move_s

function vrna_move_init(pos_5, pos_3)
    ccall((:vrna_move_init, libRNA), vrna_move_t, (Cint, Cint), pos_5, pos_3)
end

function vrna_move_list_free(moves)
    ccall((:vrna_move_list_free, libRNA), Cvoid, (Ptr{vrna_move_t},), moves)
end

function vrna_move_apply(pt, m)
    ccall((:vrna_move_apply, libRNA), Cvoid, (Ptr{Cshort}, Ptr{vrna_move_t}), pt, m)
end

function vrna_move_apply_db(structure, pt, m)
    ccall((:vrna_move_apply_db, libRNA), Cvoid, (Ptr{Cchar}, Ptr{Cshort}, Ptr{vrna_move_t}), structure, pt, m)
end

function vrna_move_is_removal(m)
    ccall((:vrna_move_is_removal, libRNA), Cint, (Ptr{vrna_move_t},), m)
end

function vrna_move_is_insertion(m)
    ccall((:vrna_move_is_insertion, libRNA), Cint, (Ptr{vrna_move_t},), m)
end

function vrna_move_is_shift(m)
    ccall((:vrna_move_is_shift, libRNA), Cint, (Ptr{vrna_move_t},), m)
end

function vrna_move_compare(a, b, pt)
    ccall((:vrna_move_compare, libRNA), Cint, (Ptr{vrna_move_t}, Ptr{vrna_move_t}, Ptr{Cshort}), a, b, pt)
end

function vrna_eval_structure(vc, structure)
    ccall((:vrna_eval_structure, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), vc, structure)
end

function vrna_eval_covar_structure(vc, structure)
    ccall((:vrna_eval_covar_structure, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), vc, structure)
end

function vrna_eval_structure_verbose(vc, structure, file)
    ccall((:vrna_eval_structure_verbose, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Ptr{Libc.FILE}), vc, structure, file)
end

function vrna_eval_structure_v(vc, structure, verbosity_level, file)
    ccall((:vrna_eval_structure_v, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), vc, structure, verbosity_level, file)
end

function vrna_eval_structure_cstr(vc, structure, verbosity_level, output_stream)
    ccall((:vrna_eval_structure_cstr, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Cint, vrna_cstr_t), vc, structure, verbosity_level, output_stream)
end

function vrna_eval_structure_pt(vc, pt)
    ccall((:vrna_eval_structure_pt, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cshort}), vc, pt)
end

function vrna_eval_structure_pt_verbose(vc, pt, file)
    ccall((:vrna_eval_structure_pt_verbose, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cshort}, Ptr{Libc.FILE}), vc, pt, file)
end

function vrna_eval_structure_pt_v(vc, pt, verbosity_level, file)
    ccall((:vrna_eval_structure_pt_v, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cshort}, Cint, Ptr{Libc.FILE}), vc, pt, verbosity_level, file)
end

function vrna_eval_structure_simple(string, structure)
    ccall((:vrna_eval_structure_simple, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), string, structure)
end

function vrna_eval_circ_structure(string, structure)
    ccall((:vrna_eval_circ_structure, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), string, structure)
end

function vrna_eval_gquad_structure(string, structure)
    ccall((:vrna_eval_gquad_structure, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), string, structure)
end

function vrna_eval_circ_gquad_structure(string, structure)
    ccall((:vrna_eval_circ_gquad_structure, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), string, structure)
end

function vrna_eval_structure_simple_verbose(string, structure, file)
    ccall((:vrna_eval_structure_simple_verbose, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Libc.FILE}), string, structure, file)
end

function vrna_eval_structure_simple_v(string, structure, verbosity_level, file)
    ccall((:vrna_eval_structure_simple_v, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), string, structure, verbosity_level, file)
end

function vrna_eval_circ_structure_v(string, structure, verbosity_level, file)
    ccall((:vrna_eval_circ_structure_v, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), string, structure, verbosity_level, file)
end

function vrna_eval_gquad_structure_v(string, structure, verbosity_level, file)
    ccall((:vrna_eval_gquad_structure_v, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), string, structure, verbosity_level, file)
end

function vrna_eval_circ_gquad_structure_v(string, structure, verbosity_level, file)
    ccall((:vrna_eval_circ_gquad_structure_v, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), string, structure, verbosity_level, file)
end

function vrna_eval_consensus_structure_simple(alignment, structure)
    ccall((:vrna_eval_consensus_structure_simple, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}), alignment, structure)
end

function vrna_eval_circ_consensus_structure(alignment, structure)
    ccall((:vrna_eval_circ_consensus_structure, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}), alignment, structure)
end

function vrna_eval_gquad_consensus_structure(alignment, structure)
    ccall((:vrna_eval_gquad_consensus_structure, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}), alignment, structure)
end

function vrna_eval_circ_gquad_consensus_structure(alignment, structure)
    ccall((:vrna_eval_circ_gquad_consensus_structure, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}), alignment, structure)
end

function vrna_eval_consensus_structure_simple_verbose(alignment, structure, file)
    ccall((:vrna_eval_consensus_structure_simple_verbose, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}, Ptr{Libc.FILE}), alignment, structure, file)
end

function vrna_eval_consensus_structure_simple_v(alignment, structure, verbosity_level, file)
    ccall((:vrna_eval_consensus_structure_simple_v, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), alignment, structure, verbosity_level, file)
end

function vrna_eval_circ_consensus_structure_v(alignment, structure, verbosity_level, file)
    ccall((:vrna_eval_circ_consensus_structure_v, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), alignment, structure, verbosity_level, file)
end

function vrna_eval_gquad_consensus_structure_v(alignment, structure, verbosity_level, file)
    ccall((:vrna_eval_gquad_consensus_structure_v, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), alignment, structure, verbosity_level, file)
end

function vrna_eval_circ_gquad_consensus_structure_v(alignment, structure, verbosity_level, file)
    ccall((:vrna_eval_circ_gquad_consensus_structure_v, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), alignment, structure, verbosity_level, file)
end

function vrna_eval_structure_pt_simple(string, pt)
    ccall((:vrna_eval_structure_pt_simple, libRNA), Cint, (Ptr{Cchar}, Ptr{Cshort}), string, pt)
end

function vrna_eval_structure_pt_simple_verbose(string, pt, file)
    ccall((:vrna_eval_structure_pt_simple_verbose, libRNA), Cint, (Ptr{Cchar}, Ptr{Cshort}, Ptr{Libc.FILE}), string, pt, file)
end

function vrna_eval_structure_pt_simple_v(string, pt, verbosity_level, file)
    ccall((:vrna_eval_structure_pt_simple_v, libRNA), Cint, (Ptr{Cchar}, Ptr{Cshort}, Cint, Ptr{Libc.FILE}), string, pt, verbosity_level, file)
end

function vrna_eval_consensus_structure_pt_simple(alignment, pt)
    ccall((:vrna_eval_consensus_structure_pt_simple, libRNA), Cint, (Ptr{Ptr{Cchar}}, Ptr{Cshort}), alignment, pt)
end

function vrna_eval_consensus_structure_pt_simple_verbose(alignment, pt, file)
    ccall((:vrna_eval_consensus_structure_pt_simple_verbose, libRNA), Cint, (Ptr{Ptr{Cchar}}, Ptr{Cshort}, Ptr{Libc.FILE}), alignment, pt, file)
end

function vrna_eval_consensus_structure_pt_simple_v(alignment, pt, verbosity_level, file)
    ccall((:vrna_eval_consensus_structure_pt_simple_v, libRNA), Cint, (Ptr{Ptr{Cchar}}, Ptr{Cshort}, Cint, Ptr{Libc.FILE}), alignment, pt, verbosity_level, file)
end

function vrna_eval_loop_pt(vc, i, pt)
    ccall((:vrna_eval_loop_pt, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cint, Ptr{Cshort}), vc, i, pt)
end

function vrna_eval_loop_pt_v(vc, i, pt, verbosity_level)
    ccall((:vrna_eval_loop_pt_v, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cint, Ptr{Cshort}, Cint), vc, i, pt, verbosity_level)
end

function vrna_eval_move(vc, structure, m1, m2)
    ccall((:vrna_eval_move, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Cint, Cint), vc, structure, m1, m2)
end

function vrna_eval_move_pt(vc, pt, m1, m2)
    ccall((:vrna_eval_move_pt, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cshort}, Cint, Cint), vc, pt, m1, m2)
end

function vrna_eval_move_pt_simple(string, pt, m1, m2)
    ccall((:vrna_eval_move_pt_simple, libRNA), Cint, (Ptr{Cchar}, Ptr{Cshort}, Cint, Cint), string, pt, m1, m2)
end

function vrna_eval_move_shift_pt(vc, m, structure)
    ccall((:vrna_eval_move_shift_pt, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{vrna_move_t}, Ptr{Cshort}), vc, m, structure)
end

function energy_of_structure(string, structure, verbosity_level)
    ccall((:energy_of_structure, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Cint), string, structure, verbosity_level)
end

function energy_of_struct_par(string, structure, parameters, verbosity_level)
    ccall((:energy_of_struct_par, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Ptr{vrna_param_t}, Cint), string, structure, parameters, verbosity_level)
end

function energy_of_circ_structure(string, structure, verbosity_level)
    ccall((:energy_of_circ_structure, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Cint), string, structure, verbosity_level)
end

function energy_of_circ_struct_par(string, structure, parameters, verbosity_level)
    ccall((:energy_of_circ_struct_par, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Ptr{vrna_param_t}, Cint), string, structure, parameters, verbosity_level)
end

function energy_of_gquad_structure(string, structure, verbosity_level)
    ccall((:energy_of_gquad_structure, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Cint), string, structure, verbosity_level)
end

function energy_of_gquad_struct_par(string, structure, parameters, verbosity_level)
    ccall((:energy_of_gquad_struct_par, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Ptr{vrna_param_t}, Cint), string, structure, parameters, verbosity_level)
end

function energy_of_structure_pt(string, ptable, s, s1, verbosity_level)
    ccall((:energy_of_structure_pt, libRNA), Cint, (Ptr{Cchar}, Ptr{Cshort}, Ptr{Cshort}, Ptr{Cshort}, Cint), string, ptable, s, s1, verbosity_level)
end

function energy_of_struct_pt_par(string, ptable, s, s1, parameters, verbosity_level)
    ccall((:energy_of_struct_pt_par, libRNA), Cint, (Ptr{Cchar}, Ptr{Cshort}, Ptr{Cshort}, Ptr{Cshort}, Ptr{vrna_param_t}, Cint), string, ptable, s, s1, parameters, verbosity_level)
end

function energy_of_move(string, structure, m1, m2)
    ccall((:energy_of_move, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Cint, Cint), string, structure, m1, m2)
end

function energy_of_move_pt(pt, s, s1, m1, m2)
    ccall((:energy_of_move_pt, libRNA), Cint, (Ptr{Cshort}, Ptr{Cshort}, Ptr{Cshort}, Cint, Cint), pt, s, s1, m1, m2)
end

function loop_energy(ptable, s, s1, i)
    ccall((:loop_energy, libRNA), Cint, (Ptr{Cshort}, Ptr{Cshort}, Ptr{Cshort}, Cint), ptable, s, s1, i)
end

function energy_of_struct(string, structure)
    ccall((:energy_of_struct, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), string, structure)
end

function energy_of_struct_pt(string, ptable, s, s1)
    ccall((:energy_of_struct_pt, libRNA), Cint, (Ptr{Cchar}, Ptr{Cshort}, Ptr{Cshort}, Ptr{Cshort}), string, ptable, s, s1)
end

function energy_of_circ_struct(string, structure)
    ccall((:energy_of_circ_struct, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), string, structure)
end

# typedef void ( vrna_callback_free_auxdata ) ( void * data )
const vrna_callback_free_auxdata = Cvoid

# typedef void ( vrna_callback_recursion_status ) ( unsigned char status , void * data )
const vrna_callback_recursion_status = Cvoid

function vrna_fold_compound(sequence, md_p, options)
    ccall((:vrna_fold_compound, libRNA), Ptr{vrna_fold_compound_t}, (Ptr{Cchar}, Ptr{vrna_md_t}, Cuint), sequence, md_p, options)
end

function vrna_fold_compound_comparative(sequences, md_p, options)
    ccall((:vrna_fold_compound_comparative, libRNA), Ptr{vrna_fold_compound_t}, (Ptr{Ptr{Cchar}}, Ptr{vrna_md_t}, Cuint), sequences, md_p, options)
end

function vrna_fold_compound_comparative2(sequences, names, orientation, start, genome_size, md_p, options)
    ccall((:vrna_fold_compound_comparative2, libRNA), Ptr{vrna_fold_compound_t}, (Ptr{Ptr{Cchar}}, Ptr{Ptr{Cchar}}, Ptr{Cuchar}, Ptr{Culonglong}, Ptr{Culonglong}, Ptr{vrna_md_t}, Cuint), sequences, names, orientation, start, genome_size, md_p, options)
end

function vrna_fold_compound_TwoD(sequence, s1, s2, md_p, options)
    ccall((:vrna_fold_compound_TwoD, libRNA), Ptr{vrna_fold_compound_t}, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Ptr{vrna_md_t}, Cuint), sequence, s1, s2, md_p, options)
end

function vrna_fold_compound_prepare(fc, options)
    ccall((:vrna_fold_compound_prepare, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cuint), fc, options)
end

function vrna_fold_compound_free(fc)
    ccall((:vrna_fold_compound_free, libRNA), Cvoid, (Ptr{vrna_fold_compound_t},), fc)
end

function vrna_fold_compound_add_auxdata(fc, data, f)
    ccall((:vrna_fold_compound_add_auxdata, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}, Ptr{Cvoid}), fc, data, f)
end

function vrna_fold_compound_add_callback(fc, f)
    ccall((:vrna_fold_compound_add_callback, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Ptr{Cvoid}), fc, f)
end

# typedef void ( vrna_heat_capacity_callback ) ( float temp , float heat_capacity , void * data )
const vrna_heat_capacity_callback = Cvoid

struct vrna_heat_capacity_s
    temperature::Cfloat
    heat_capacity::Cfloat
end

const vrna_heat_capacity_t = vrna_heat_capacity_s

function vrna_heat_capacity(fc, T_min, T_max, T_increment, mpoints)
    ccall((:vrna_heat_capacity, libRNA), Ptr{vrna_heat_capacity_t}, (Ptr{vrna_fold_compound_t}, Cfloat, Cfloat, Cfloat, Cuint), fc, T_min, T_max, T_increment, mpoints)
end

function vrna_heat_capacity_cb(fc, T_min, T_max, T_increment, mpoints, cb, data)
    ccall((:vrna_heat_capacity_cb, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Cfloat, Cfloat, Cfloat, Cuint, Ptr{Cvoid}, Ptr{Cvoid}), fc, T_min, T_max, T_increment, mpoints, cb, data)
end

function vrna_heat_capacity_simple(sequence, T_min, T_max, T_increment, mpoints)
    ccall((:vrna_heat_capacity_simple, libRNA), Ptr{vrna_heat_capacity_t}, (Ptr{Cchar}, Cfloat, Cfloat, Cfloat, Cuint), sequence, T_min, T_max, T_increment, mpoints)
end

function inverse_fold(start, target)
    ccall((:inverse_fold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), start, target)
end

function inverse_pf_fold(start, target)
    ccall((:inverse_pf_fold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), start, target)
end

# typedef void ( vrna_callback_move_update ) ( vrna_fold_compound_t * fc , vrna_move_t neighbor , unsigned int state , void * data )
const vrna_callback_move_update = Cvoid

function vrna_loopidx_update(loopidx, pt, length, m)
    ccall((:vrna_loopidx_update, libRNA), Cvoid, (Ptr{Cint}, Ptr{Cshort}, Cint, Ptr{vrna_move_t}), loopidx, pt, length, m)
end

function vrna_neighbors(vc, pt, options)
    ccall((:vrna_neighbors, libRNA), Ptr{vrna_move_t}, (Ptr{vrna_fold_compound_t}, Ptr{Cshort}, Cuint), vc, pt, options)
end

function vrna_neighbors_successive(vc, curr_move, prev_pt, prev_neighbors, size_prev_neighbors, size_neighbors, options)
    ccall((:vrna_neighbors_successive, libRNA), Ptr{vrna_move_t}, (Ptr{vrna_fold_compound_t}, Ptr{vrna_move_t}, Ptr{Cshort}, Ptr{vrna_move_t}, Cint, Ptr{Cint}, Cuint), vc, curr_move, prev_pt, prev_neighbors, size_prev_neighbors, size_neighbors, options)
end

function vrna_move_neighbor_diff_cb(fc, ptable, move, cb, data, options)
    ccall((:vrna_move_neighbor_diff_cb, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cshort}, vrna_move_t, Ptr{Cvoid}, Ptr{Cvoid}, Cuint), fc, ptable, move, cb, data, options)
end

function vrna_move_neighbor_diff(fc, ptable, move, invalid_moves, options)
    ccall((:vrna_move_neighbor_diff, libRNA), Ptr{vrna_move_t}, (Ptr{vrna_fold_compound_t}, Ptr{Cshort}, vrna_move_t, Ptr{Ptr{vrna_move_t}}, Cuint), fc, ptable, move, invalid_moves, options)
end

function vrna_MEA(fc, gamma, mea)
    ccall((:vrna_MEA, libRNA), Ptr{Cchar}, (Ptr{vrna_fold_compound_t}, Cdouble, Ptr{Cfloat}), fc, gamma, mea)
end

function vrna_MEA_from_plist(plist_, sequence, gamma, md, mea)
    ccall((:vrna_MEA_from_plist, libRNA), Ptr{Cchar}, (Ptr{vrna_ep_t}, Ptr{Cchar}, Cdouble, Ptr{vrna_md_t}, Ptr{Cfloat}), plist_, sequence, gamma, md, mea)
end

function MEA(p, structure, gamma)
    ccall((:MEA, libRNA), Cfloat, (Ptr{plist}, Ptr{Cchar}, Cdouble), p, structure, gamma)
end

function MEA_seq(p, sequence, structure, gamma, pf)
    ccall((:MEA_seq, libRNA), Cfloat, (Ptr{plist}, Ptr{Cchar}, Ptr{Cchar}, Cdouble, Ptr{vrna_exp_param_t}), p, sequence, structure, gamma, pf)
end

function vrna_mfe(vc, structure)
    ccall((:vrna_mfe, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), vc, structure)
end

function vrna_mfe_dimer(vc, structure)
    ccall((:vrna_mfe_dimer, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), vc, structure)
end

function vrna_fold(sequence, structure)
    ccall((:vrna_fold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), sequence, structure)
end

function vrna_circfold(sequence, structure)
    ccall((:vrna_circfold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), sequence, structure)
end

function vrna_alifold(sequences, structure)
    ccall((:vrna_alifold, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}), sequences, structure)
end

function vrna_circalifold(sequences, structure)
    ccall((:vrna_circalifold, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}), sequences, structure)
end

function vrna_cofold(sequence, structure)
    ccall((:vrna_cofold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), sequence, structure)
end

function vrna_backtrack_from_intervals(vc, bp_stack, bt_stack, s)
    ccall((:vrna_backtrack_from_intervals, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{vrna_bp_stack_t}, Ptr{sect}, Cint), vc, bp_stack, bt_stack, s)
end

function vrna_backtrack5(fc, length, structure)
    ccall((:vrna_backtrack5, libRNA), Cfloat, (Ptr{vrna_fold_compound_t}, Cuint, Ptr{Cchar}), fc, length, structure)
end

function vrna_backtrack_window(fc, Lfold_filename, file_pos, structure, mfe)
    ccall((:vrna_backtrack_window, libRNA), Cint, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}, Clong, Ptr{Ptr{Cchar}}, Cdouble), fc, Lfold_filename, file_pos, structure, mfe)
end

function vrna_params_load(fname, options)
    ccall((:vrna_params_load, libRNA), Cint, (Ptr{Cchar}, Cuint), fname, options)
end

function vrna_params_save(fname, options)
    ccall((:vrna_params_save, libRNA), Cint, (Ptr{Cchar}, Cuint), fname, options)
end

function vrna_params_load_from_string(string, name, options)
    ccall((:vrna_params_load_from_string, libRNA), Cint, (Ptr{Cchar}, Ptr{Cchar}, Cuint), string, name, options)
end

function vrna_params_load_defaults()
    ccall((:vrna_params_load_defaults, libRNA), Cint, ())
end

function vrna_params_load_RNA_Turner2004()
    ccall((:vrna_params_load_RNA_Turner2004, libRNA), Cint, ())
end

function vrna_params_load_RNA_Turner1999()
    ccall((:vrna_params_load_RNA_Turner1999, libRNA), Cint, ())
end

function vrna_params_load_RNA_Andronescu2007()
    ccall((:vrna_params_load_RNA_Andronescu2007, libRNA), Cint, ())
end

function vrna_params_load_RNA_Langdon2018()
    ccall((:vrna_params_load_RNA_Langdon2018, libRNA), Cint, ())
end

function vrna_params_load_RNA_misc_special_hairpins()
    ccall((:vrna_params_load_RNA_misc_special_hairpins, libRNA), Cint, ())
end

function vrna_params_load_DNA_Mathews2004()
    ccall((:vrna_params_load_DNA_Mathews2004, libRNA), Cint, ())
end

function vrna_params_load_DNA_Mathews1999()
    ccall((:vrna_params_load_DNA_Mathews1999, libRNA), Cint, ())
end

@cenum parset::Int32 begin
    UNKNOWN = -1
    QUIT = 0
    S = 1
    S_H = 2
    HP = 3
    HP_H = 4
    B = 5
    B_H = 6
    IL = 7
    IL_H = 8
    MMH = 9
    MMH_H = 10
    MMI = 11
    MMI_H = 12
    MMI1N = 13
    MMI1N_H = 14
    MMI23 = 15
    MMI23_H = 16
    MMM = 17
    MMM_H = 18
    MME = 19
    MME_H = 20
    D5 = 21
    D5_H = 22
    D3 = 23
    D3_H = 24
    INT11 = 25
    INT11_H = 26
    INT21 = 27
    INT21_H = 28
    INT22 = 29
    INT22_H = 30
    ML = 31
    TL = 32
    TRI = 33
    HEX = 34
    NIN = 35
    MISC = 36
end

function last_parameter_file()
    ccall((:last_parameter_file, libRNA), Ptr{Cchar}, ())
end

function read_parameter_file(fname)
    ccall((:read_parameter_file, libRNA), Cvoid, (Ptr{Cchar},), fname)
end

function write_parameter_file(fname)
    ccall((:write_parameter_file, libRNA), Cvoid, (Ptr{Cchar},), fname)
end

function gettype(ident)
    ccall((:gettype, libRNA), parset, (Ptr{Cchar},), ident)
end

function settype(s)
    ccall((:settype, libRNA), Ptr{Cchar}, (parset,), s)
end

struct vrna_dimer_pf_s
    F0AB::Cdouble
    FAB::Cdouble
    FcAB::Cdouble
    FA::Cdouble
    FB::Cdouble
end

const vrna_dimer_pf_t = vrna_dimer_pf_s

struct vrna_multimer_pf_s
    F_connected::Cdouble
    F_monomers::Ptr{Cdouble}
    num_monomers::Csize_t
end

const vrna_multimer_pf_t = vrna_multimer_pf_s

const cofoldF = vrna_dimer_pf_s

function vrna_centroid(vc, dist)
    ccall((:vrna_centroid, libRNA), Ptr{Cchar}, (Ptr{vrna_fold_compound_t}, Ptr{Cdouble}), vc, dist)
end

function vrna_centroid_from_plist(length, dist, pl)
    ccall((:vrna_centroid_from_plist, libRNA), Ptr{Cchar}, (Cint, Ptr{Cdouble}, Ptr{vrna_ep_t}), length, dist, pl)
end

function vrna_centroid_from_probs(length, dist, probs)
    ccall((:vrna_centroid_from_probs, libRNA), Ptr{Cchar}, (Cint, Ptr{Cdouble}, Ptr{FLT_OR_DBL}), length, dist, probs)
end

function get_centroid_struct_pl(length, dist, pl)
    ccall((:get_centroid_struct_pl, libRNA), Ptr{Cchar}, (Cint, Ptr{Cdouble}, Ptr{vrna_ep_t}), length, dist, pl)
end

function get_centroid_struct_pr(length, dist, pr)
    ccall((:get_centroid_struct_pr, libRNA), Ptr{Cchar}, (Cint, Ptr{Cdouble}, Ptr{FLT_OR_DBL}), length, dist, pr)
end

# typedef void ( vrna_boltzmann_sampling_callback ) ( const char * stucture , void * data )
const vrna_boltzmann_sampling_callback = Cvoid

mutable struct vrna_pbacktrack_memory_s end

const vrna_pbacktrack_mem_t = Ptr{vrna_pbacktrack_memory_s}

function vrna_pbacktrack5(fc, length)
    ccall((:vrna_pbacktrack5, libRNA), Ptr{Cchar}, (Ptr{vrna_fold_compound_t}, Cuint), fc, length)
end

function vrna_pbacktrack5_num(fc, num_samples, length, options)
    ccall((:vrna_pbacktrack5_num, libRNA), Ptr{Ptr{Cchar}}, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Cuint), fc, num_samples, length, options)
end

function vrna_pbacktrack5_cb(fc, num_samples, length, cb, data, options)
    ccall((:vrna_pbacktrack5_cb, libRNA), Cuint, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Ptr{Cvoid}, Ptr{Cvoid}, Cuint), fc, num_samples, length, cb, data, options)
end

function vrna_pbacktrack5_resume(vc, num_samples, length, nr_mem, options)
    ccall((:vrna_pbacktrack5_resume, libRNA), Ptr{Ptr{Cchar}}, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Ptr{vrna_pbacktrack_mem_t}, Cuint), vc, num_samples, length, nr_mem, options)
end

function vrna_pbacktrack5_resume_cb(fc, num_samples, length, cb, data, nr_mem, options)
    ccall((:vrna_pbacktrack5_resume_cb, libRNA), Cuint, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{vrna_pbacktrack_mem_t}, Cuint), fc, num_samples, length, cb, data, nr_mem, options)
end

function vrna_pbacktrack(fc)
    ccall((:vrna_pbacktrack, libRNA), Ptr{Cchar}, (Ptr{vrna_fold_compound_t},), fc)
end

function vrna_pbacktrack_num(fc, num_samples, options)
    ccall((:vrna_pbacktrack_num, libRNA), Ptr{Ptr{Cchar}}, (Ptr{vrna_fold_compound_t}, Cuint, Cuint), fc, num_samples, options)
end

function vrna_pbacktrack_cb(fc, num_samples, cb, data, options)
    ccall((:vrna_pbacktrack_cb, libRNA), Cuint, (Ptr{vrna_fold_compound_t}, Cuint, Ptr{Cvoid}, Ptr{Cvoid}, Cuint), fc, num_samples, cb, data, options)
end

function vrna_pbacktrack_resume(fc, num_samples, nr_mem, options)
    ccall((:vrna_pbacktrack_resume, libRNA), Ptr{Ptr{Cchar}}, (Ptr{vrna_fold_compound_t}, Cuint, Ptr{vrna_pbacktrack_mem_t}, Cuint), fc, num_samples, nr_mem, options)
end

function vrna_pbacktrack_resume_cb(fc, num_samples, cb, data, nr_mem, options)
    ccall((:vrna_pbacktrack_resume_cb, libRNA), Cuint, (Ptr{vrna_fold_compound_t}, Cuint, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{vrna_pbacktrack_mem_t}, Cuint), fc, num_samples, cb, data, nr_mem, options)
end

function vrna_pbacktrack_sub(fc, start, _end)
    ccall((:vrna_pbacktrack_sub, libRNA), Ptr{Cchar}, (Ptr{vrna_fold_compound_t}, Cuint, Cuint), fc, start, _end)
end

function vrna_pbacktrack_sub_num(fc, num_samples, start, _end, options)
    ccall((:vrna_pbacktrack_sub_num, libRNA), Ptr{Ptr{Cchar}}, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Cuint, Cuint), fc, num_samples, start, _end, options)
end

function vrna_pbacktrack_sub_cb(fc, num_samples, start, _end, cb, data, options)
    ccall((:vrna_pbacktrack_sub_cb, libRNA), Cuint, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Cuint, Ptr{Cvoid}, Ptr{Cvoid}, Cuint), fc, num_samples, start, _end, cb, data, options)
end

function vrna_pbacktrack_sub_resume(vc, num_samples, start, _end, nr_mem, options)
    ccall((:vrna_pbacktrack_sub_resume, libRNA), Ptr{Ptr{Cchar}}, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Cuint, Ptr{vrna_pbacktrack_mem_t}, Cuint), vc, num_samples, start, _end, nr_mem, options)
end

function vrna_pbacktrack_sub_resume_cb(fc, num_samples, start, _end, bs_cb, data, nr_mem, options)
    ccall((:vrna_pbacktrack_sub_resume_cb, libRNA), Cuint, (Ptr{vrna_fold_compound_t}, Cuint, Cuint, Cuint, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{vrna_pbacktrack_mem_t}, Cuint), fc, num_samples, start, _end, bs_cb, data, nr_mem, options)
end

function vrna_pbacktrack_mem_free(s)
    ccall((:vrna_pbacktrack_mem_free, libRNA), Cvoid, (vrna_pbacktrack_mem_t,), s)
end

function vrna_pf(vc, structure)
    ccall((:vrna_pf, libRNA), FLT_OR_DBL, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), vc, structure)
end

function vrna_pf_dimer(vc, structure)
    ccall((:vrna_pf_dimer, libRNA), vrna_dimer_pf_t, (Ptr{vrna_fold_compound_t}, Ptr{Cchar}), vc, structure)
end

function vrna_pf_substrands(fc, complex_size)
    ccall((:vrna_pf_substrands, libRNA), Ptr{FLT_OR_DBL}, (Ptr{vrna_fold_compound_t}, Csize_t), fc, complex_size)
end

function vrna_pf_add(dG1, dG2, kT)
    ccall((:vrna_pf_add, libRNA), FLT_OR_DBL, (FLT_OR_DBL, FLT_OR_DBL, Cdouble), dG1, dG2, kT)
end

function vrna_pf_fold(sequence, structure, pl)
    ccall((:vrna_pf_fold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{vrna_ep_t}}), sequence, structure, pl)
end

function vrna_pf_circfold(sequence, structure, pl)
    ccall((:vrna_pf_circfold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{vrna_ep_t}}), sequence, structure, pl)
end

function vrna_pf_alifold(sequences, structure, pl)
    ccall((:vrna_pf_alifold, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}, Ptr{Ptr{vrna_ep_t}}), sequences, structure, pl)
end

function vrna_pf_circalifold(sequences, structure, pl)
    ccall((:vrna_pf_circalifold, libRNA), Cfloat, (Ptr{Ptr{Cchar}}, Ptr{Cchar}, Ptr{Ptr{vrna_ep_t}}), sequences, structure, pl)
end

function vrna_pf_co_fold(seq, structure, pl)
    ccall((:vrna_pf_co_fold, libRNA), vrna_dimer_pf_t, (Ptr{Cchar}, Ptr{Cchar}, Ptr{Ptr{vrna_ep_t}}), seq, structure, pl)
end

function vrna_pf_float_precision()
    ccall((:vrna_pf_float_precision, libRNA), Cint, ())
end

function pf_fold_par(sequence, structure, parameters, calculate_bppm, is_constrained, is_circular)
    ccall((:pf_fold_par, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}, Ptr{vrna_exp_param_t}, Cint, Cint, Cint), sequence, structure, parameters, calculate_bppm, is_constrained, is_circular)
end

function pf_fold(sequence, structure)
    ccall((:pf_fold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), sequence, structure)
end

function pf_circ_fold(sequence, structure)
    ccall((:pf_circ_fold, libRNA), Cfloat, (Ptr{Cchar}, Ptr{Cchar}), sequence, structure)
end

function pbacktrack(sequence)
    ccall((:pbacktrack, libRNA), Ptr{Cchar}, (Ptr{Cchar},), sequence)
end

function pbacktrack5(sequence, length)
    ccall((:pbacktrack5, libRNA), Ptr{Cchar}, (Ptr{Cchar}, Cint), sequence, length)
end

function pbacktrack_circ(sequence)
    ccall((:pbacktrack_circ, libRNA), Ptr{Cchar}, (Ptr{Cchar},), sequence)
end

function free_pf_arrays()
    ccall((:free_pf_arrays, libRNA), Cvoid, ())
end

function update_pf_params(length)
    ccall((:update_pf_params, libRNA), Cvoid, (Cint,), length)
end

function update_pf_params_par(length, parameters)
    ccall((:update_pf_params_par, libRNA), Cvoid, (Cint, Ptr{vrna_exp_param_t}), length, parameters)
end

function export_bppm()
    ccall((:export_bppm, libRNA), Ptr{FLT_OR_DBL}, ())
end

function get_pf_arrays(S_p, S1_p, ptype_p, qb_p, qm_p, q1k_p, qln_p)
    ccall((:get_pf_arrays, libRNA), Cint, (Ptr{Ptr{Cshort}}, Ptr{Ptr{Cshort}}, Ptr{Ptr{Cchar}}, Ptr{Ptr{FLT_OR_DBL}}, Ptr{Ptr{FLT_OR_DBL}}, Ptr{Ptr{FLT_OR_DBL}}, Ptr{Ptr{FLT_OR_DBL}}), S_p, S1_p, ptype_p, qb_p, qm_p, q1k_p, qln_p)
end

function get_subseq_F(i, j)
    ccall((:get_subseq_F, libRNA), Cdouble, (Cint, Cint), i, j)
end

function mean_bp_distance(length)
    ccall((:mean_bp_distance, libRNA), Cdouble, (Cint,), length)
end

function mean_bp_distance_pr(length, pr)
    ccall((:mean_bp_distance_pr, libRNA), Cdouble, (Cint, Ptr{FLT_OR_DBL}), length, pr)
end

function stackProb(cutoff)
    ccall((:stackProb, libRNA), Ptr{vrna_ep_t}, (Cdouble,), cutoff)
end

function init_pf_fold(length)
    ccall((:init_pf_fold, libRNA), Cvoid, (Cint,), length)
end

function centroid(length, dist)
    ccall((:centroid, libRNA), Ptr{Cchar}, (Cint, Ptr{Cdouble}), length, dist)
end

function get_centroid_struct_gquad_pr(length, dist)
    ccall((:get_centroid_struct_gquad_pr, libRNA), Ptr{Cchar}, (Cint, Ptr{Cdouble}), length, dist)
end

function mean_bp_dist(length)
    ccall((:mean_bp_dist, libRNA), Cdouble, (Cint,), length)
end

function expLoopEnergy(u1, u2, type, type2, si1, sj1, sp1, sq1)
    ccall((:expLoopEnergy, libRNA), Cdouble, (Cint, Cint, Cint, Cint, Cshort, Cshort, Cshort, Cshort), u1, u2, type, type2, si1, sj1, sp1, sq1)
end

function expHairpinEnergy(u, type, si1, sj1, string)
    ccall((:expHairpinEnergy, libRNA), Cdouble, (Cint, Cint, Cshort, Cshort, Ptr{Cchar}), u, type, si1, sj1, string)
end

function assign_plist_gquad_from_pr(pl, length, cut_off)
    ccall((:assign_plist_gquad_from_pr, libRNA), Cvoid, (Ptr{Ptr{vrna_ep_t}}, Cint, Cdouble), pl, length, cut_off)
end

struct vrna_plot_layout_s
    length::Cuint
    x::Ptr{Cfloat}
    y::Ptr{Cfloat}
    arcs::Ptr{Cdouble}
    bbox::NTuple{4, Cint}
end

const vrna_plot_layout_t = vrna_plot_layout_s

function vrna_plot_coords_naview(structure, x, y)
    ccall((:vrna_plot_coords_naview, libRNA), Cint, (Ptr{Cchar}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}), structure, x, y)
end

function vrna_plot_coords_naview_pt(pt, x, y)
    ccall((:vrna_plot_coords_naview_pt, libRNA), Cint, (Ptr{Cshort}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}), pt, x, y)
end

function naview_xy_coordinates(pair_table, X, Y)
    ccall((:naview_xy_coordinates, libRNA), Cint, (Ptr{Cshort}, Ptr{Cfloat}, Ptr{Cfloat}), pair_table, X, Y)
end

function vrna_plot_coords_turtle(structure, x, y, arc_coords)
    ccall((:vrna_plot_coords_turtle, libRNA), Cint, (Ptr{Cchar}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cdouble}}), structure, x, y, arc_coords)
end

function vrna_plot_coords_turtle_pt(pair_table, x, y, arc_coords)
    ccall((:vrna_plot_coords_turtle_pt, libRNA), Cint, (Ptr{Cshort}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cdouble}}), pair_table, x, y, arc_coords)
end

struct vrna_plot_options_puzzler_t
    drawArcs::Cshort
    paired::Cdouble
    unpaired::Cdouble
    checkAncestorIntersections::Cshort
    checkSiblingIntersections::Cshort
    checkExteriorIntersections::Cshort
    allowFlipping::Cshort
    optimize::Cshort
    maximumNumberOfConfigChangesAllowed::Cint
    config::Ptr{Cchar}
    filename::Ptr{Cchar}
    numberOfChangesAppliedToConfig::Cint
    psNumber::Cint
end

function vrna_plot_coords_puzzler(structure, x, y, arc_coords, options)
    ccall((:vrna_plot_coords_puzzler, libRNA), Cint, (Ptr{Cchar}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cdouble}}, Ptr{vrna_plot_options_puzzler_t}), structure, x, y, arc_coords, options)
end

function vrna_plot_coords_puzzler_pt(pair_table, x, y, arc_coords, puzzler)
    ccall((:vrna_plot_coords_puzzler_pt, libRNA), Cint, (Ptr{Cshort}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cdouble}}, Ptr{vrna_plot_options_puzzler_t}), pair_table, x, y, arc_coords, puzzler)
end

function vrna_plot_options_puzzler()
    ccall((:vrna_plot_options_puzzler, libRNA), Ptr{vrna_plot_options_puzzler_t}, ())
end

function vrna_plot_options_puzzler_free(options)
    ccall((:vrna_plot_options_puzzler_free, libRNA), Cvoid, (Ptr{vrna_plot_options_puzzler_t},), options)
end

function vrna_plot_layout(structure, plot_type)
    ccall((:vrna_plot_layout, libRNA), Ptr{vrna_plot_layout_t}, (Ptr{Cchar}, Cuint), structure, plot_type)
end

function vrna_plot_layout_simple(structure)
    ccall((:vrna_plot_layout_simple, libRNA), Ptr{vrna_plot_layout_t}, (Ptr{Cchar},), structure)
end

function vrna_plot_layout_naview(structure)
    ccall((:vrna_plot_layout_naview, libRNA), Ptr{vrna_plot_layout_t}, (Ptr{Cchar},), structure)
end

function vrna_plot_layout_circular(structure)
    ccall((:vrna_plot_layout_circular, libRNA), Ptr{vrna_plot_layout_t}, (Ptr{Cchar},), structure)
end

function vrna_plot_layout_turtle(structure)
    ccall((:vrna_plot_layout_turtle, libRNA), Ptr{vrna_plot_layout_t}, (Ptr{Cchar},), structure)
end

function vrna_plot_layout_puzzler(structure, options)
    ccall((:vrna_plot_layout_puzzler, libRNA), Ptr{vrna_plot_layout_t}, (Ptr{Cchar}, Ptr{vrna_plot_options_puzzler_t}), structure, options)
end

function vrna_plot_layout_free(layout)
    ccall((:vrna_plot_layout_free, libRNA), Cvoid, (Ptr{vrna_plot_layout_t},), layout)
end

function vrna_plot_coords(structure, x, y, plot_type)
    ccall((:vrna_plot_coords, libRNA), Cint, (Ptr{Cchar}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}, Cint), structure, x, y, plot_type)
end

function vrna_plot_coords_pt(pt, x, y, plot_type)
    ccall((:vrna_plot_coords_pt, libRNA), Cint, (Ptr{Cshort}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}, Cint), pt, x, y, plot_type)
end

function vrna_plot_coords_simple(structure, x, y)
    ccall((:vrna_plot_coords_simple, libRNA), Cint, (Ptr{Cchar}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}), structure, x, y)
end

function vrna_plot_coords_simple_pt(pt, x, y)
    ccall((:vrna_plot_coords_simple_pt, libRNA), Cint, (Ptr{Cshort}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}), pt, x, y)
end

function vrna_plot_coords_circular(structure, x, y)
    ccall((:vrna_plot_coords_circular, libRNA), Cint, (Ptr{Cchar}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}), structure, x, y)
end

function vrna_plot_coords_circular_pt(pt, x, y)
    ccall((:vrna_plot_coords_circular_pt, libRNA), Cint, (Ptr{Cshort}, Ptr{Ptr{Cfloat}}, Ptr{Ptr{Cfloat}}), pt, x, y)
end

struct COORDINATE
    X::Cfloat
    Y::Cfloat
end

function simple_xy_coordinates(pair_table, X, Y)
    ccall((:simple_xy_coordinates, libRNA), Cint, (Ptr{Cshort}, Ptr{Cfloat}, Ptr{Cfloat}), pair_table, X, Y)
end

function simple_circplot_coordinates(pair_table, x, y)
    ccall((:simple_circplot_coordinates, libRNA), Cint, (Ptr{Cshort}, Ptr{Cfloat}, Ptr{Cfloat}), pair_table, x, y)
end

struct vrna_subopt_sol_s
    energy::Cfloat
    structure::Ptr{Cchar}
end

const vrna_subopt_solution_t = vrna_subopt_sol_s

# typedef void ( vrna_subopt_callback ) ( const char * stucture , float energy , void * data )
const vrna_subopt_callback = Cvoid

const SOLUTION = vrna_subopt_sol_s

function vrna_subopt(fc, delta, sorted, fp)
    ccall((:vrna_subopt, libRNA), Ptr{vrna_subopt_solution_t}, (Ptr{vrna_fold_compound_t}, Cint, Cint, Ptr{Libc.FILE}), fc, delta, sorted, fp)
end

function vrna_subopt_cb(fc, delta, cb, data)
    ccall((:vrna_subopt_cb, libRNA), Cvoid, (Ptr{vrna_fold_compound_t}, Cint, Ptr{Cvoid}, Ptr{Cvoid}), fc, delta, cb, data)
end

function subopt(seq, structure, delta, fp)
    ccall((:subopt, libRNA), Ptr{SOLUTION}, (Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), seq, structure, delta, fp)
end

function subopt_par(seq, structure, parameters, delta, is_constrained, is_circular, fp)
    ccall((:subopt_par, libRNA), Ptr{SOLUTION}, (Ptr{Cchar}, Ptr{Cchar}, Ptr{vrna_param_t}, Cint, Cint, Cint, Ptr{Libc.FILE}), seq, structure, parameters, delta, is_constrained, is_circular, fp)
end

function subopt_circ(seq, sequence, delta, fp)
    ccall((:subopt_circ, libRNA), Ptr{SOLUTION}, (Ptr{Cchar}, Ptr{Cchar}, Cint, Ptr{Libc.FILE}), seq, sequence, delta, fp)
end

function zukersubopt(string)
    ccall((:zukersubopt, libRNA), Ptr{SOLUTION}, (Ptr{Cchar},), string)
end

function zukersubopt_par(string, parameters)
    ccall((:zukersubopt_par, libRNA), Ptr{SOLUTION}, (Ptr{Cchar}, Ptr{vrna_param_t}), string, parameters)
end

function vrna_subopt_zuker(fc)
    ccall((:vrna_subopt_zuker, libRNA), Ptr{vrna_subopt_solution_t}, (Ptr{vrna_fold_compound_t},), fc)
end

struct Postorder_list
    type::Cint
    weight::Cint
    father::Cint
    sons::Cint
    leftmostleaf::Cint
end

struct Tree
    postorder_list::Ptr{Postorder_list}
    keyroots::Ptr{Cint}
end

struct swString
    type::Cint
    sign::Cint
    weight::Cfloat
end

function make_tree(struc)
    ccall((:make_tree, libRNA), Ptr{Tree}, (Ptr{Cchar},), struc)
end

function tree_edit_distance(T1, T2)
    ccall((:tree_edit_distance, libRNA), Cfloat, (Ptr{Tree}, Ptr{Tree}), T1, T2)
end

function print_tree(t)
    ccall((:print_tree, libRNA), Cvoid, (Ptr{Tree},), t)
end

function free_tree(t)
    ccall((:free_tree, libRNA), Cvoid, (Ptr{Tree},), t)
end

function vrna_alloc(size)
    ccall((:vrna_alloc, libRNA), Ptr{Cvoid}, (Cuint,), size)
end

function vrna_realloc(p, size)
    ccall((:vrna_realloc, libRNA), Ptr{Cvoid}, (Ptr{Cvoid}, Cuint), p, size)
end

function vrna_init_rand()
    ccall((:vrna_init_rand, libRNA), Cvoid, ())
end

function vrna_urn()
    ccall((:vrna_urn, libRNA), Cdouble, ())
end

function vrna_int_urn(from, to)
    ccall((:vrna_int_urn, libRNA), Cint, (Cint, Cint), from, to)
end

function vrna_time_stamp()
    ccall((:vrna_time_stamp, libRNA), Ptr{Cchar}, ())
end

function get_input_line(string, options)
    ccall((:get_input_line, libRNA), Cuint, (Ptr{Ptr{Cchar}}, Cuint), string, options)
end

function vrna_idx_row_wise(length)
    ccall((:vrna_idx_row_wise, libRNA), Ptr{Cint}, (Cuint,), length)
end

function vrna_idx_col_wise(length)
    ccall((:vrna_idx_col_wise, libRNA), Ptr{Cint}, (Cuint,), length)
end

function vrna_message_input_seq_simple()
    ccall((:vrna_message_input_seq_simple, libRNA), Cvoid, ())
end

function vrna_message_input_seq(s)
    ccall((:vrna_message_input_seq, libRNA), Cvoid, (Ptr{Cchar},), s)
end

function vrna_message_input_msa(s)
    ccall((:vrna_message_input_msa, libRNA), Cvoid, (Ptr{Cchar},), s)
end

function get_indx(length)
    ccall((:get_indx, libRNA), Ptr{Cint}, (Cuint,), length)
end

function get_iindx(length)
    ccall((:get_iindx, libRNA), Ptr{Cint}, (Cuint,), length)
end

function get_line(fp)
    ccall((:get_line, libRNA), Ptr{Cchar}, (Ptr{Libc.FILE},), fp)
end

function print_tty_input_seq()
    ccall((:print_tty_input_seq, libRNA), Cvoid, ())
end

function print_tty_input_seq_str(s)
    ccall((:print_tty_input_seq_str, libRNA), Cvoid, (Ptr{Cchar},), s)
end

function warn_user(message)
    ccall((:warn_user, libRNA), Cvoid, (Ptr{Cchar},), message)
end

function nrerror(message)
    ccall((:nrerror, libRNA), Cvoid, (Ptr{Cchar},), message)
end

function space(size)
    ccall((:space, libRNA), Ptr{Cvoid}, (Cuint,), size)
end

function xrealloc(p, size)
    ccall((:xrealloc, libRNA), Ptr{Cvoid}, (Ptr{Cvoid}, Cuint), p, size)
end

function init_rand()
    ccall((:init_rand, libRNA), Cvoid, ())
end

function urn()
    ccall((:urn, libRNA), Cdouble, ())
end

function int_urn(from, to)
    ccall((:int_urn, libRNA), Cint, (Cint, Cint), from, to)
end

function filecopy(from, to)
    ccall((:filecopy, libRNA), Cvoid, (Ptr{Libc.FILE}, Ptr{Libc.FILE}), from, to)
end

function time_stamp()
    ccall((:time_stamp, libRNA), Ptr{Cchar}, ())
end

const GASCONST = 1.98717

const K0 = 273.15

const INF = 10000000

const EMAX = INF ÷ 10

const FORBIDDEN = 9999

const BONUS = 10000

const NBPAIRS = 7

const TURN = 3

const MAXLOOP = 30

const UNIT = 100

const MINPSCORE = -2 * UNIT

const NBASES = 8

const VRNA_MODEL_DEFAULT_TEMPERATURE = 37.0

const VRNA_MODEL_DEFAULT_PF_SCALE = -1

const VRNA_MODEL_DEFAULT_BETA_SCALE = 1.0

const VRNA_MODEL_DEFAULT_DANGLES = 2

const VRNA_MODEL_DEFAULT_SPECIAL_HP = 1

const VRNA_MODEL_DEFAULT_NO_LP = 0

const VRNA_MODEL_DEFAULT_NO_GU = 0

const VRNA_MODEL_DEFAULT_NO_GU_CLOSURE = 0

const VRNA_MODEL_DEFAULT_CIRC = 0

const VRNA_MODEL_DEFAULT_GQUAD = 0

const VRNA_MODEL_DEFAULT_UNIQ_ML = 0

const VRNA_MODEL_DEFAULT_ENERGY_SET = 0

const VRNA_MODEL_DEFAULT_BACKTRACK = 1

const VRNA_MODEL_DEFAULT_BACKTRACK_TYPE = Cchar('F')

const VRNA_MODEL_DEFAULT_COMPUTE_BPP = 1

const VRNA_MODEL_DEFAULT_MAX_BP_SPAN = -1

const VRNA_MODEL_DEFAULT_WINDOW_SIZE = -1

const VRNA_MODEL_DEFAULT_LOG_ML = 0

const VRNA_MODEL_DEFAULT_ALI_OLD_EN = 0

const VRNA_MODEL_DEFAULT_ALI_RIBO = 0

const VRNA_MODEL_DEFAULT_ALI_CV_FACT = 1.0

const VRNA_MODEL_DEFAULT_ALI_NC_FACT = 1.0

const VRNA_MODEL_DEFAULT_PF_SMOOTH = 1

const MAXALPHA = 20

const model_detailsT = vrna_md_t

const VRNA_GQUAD_MAX_STACK_SIZE = 7

const VRNA_GQUAD_MIN_STACK_SIZE = 2

const VRNA_GQUAD_MAX_LINKER_LENGTH = 15

const VRNA_GQUAD_MIN_LINKER_LENGTH = 1

const VRNA_GQUAD_MIN_BOX_SIZE = 4VRNA_GQUAD_MIN_STACK_SIZE + 3VRNA_GQUAD_MIN_LINKER_LENGTH

const VRNA_GQUAD_MAX_BOX_SIZE = 4VRNA_GQUAD_MAX_STACK_SIZE + 3VRNA_GQUAD_MAX_LINKER_LENGTH

const VRNA_SEQUENCE_RNA = Cuint(1)

const VRNA_SEQUENCE_DNA = Cuint(2)

const VRNA_CONSTRAINT_FILE = 0

const VRNA_CONSTRAINT_SOFT_MFE = 0

const VRNA_OPTION_PF = Cuint(2)

const VRNA_CONSTRAINT_SOFT_PF = VRNA_OPTION_PF

const VRNA_DECOMP_PAIR_HP = Cuchar(1)

const VRNA_DECOMP_PAIR_IL = Cuchar(2)

const VRNA_DECOMP_PAIR_ML = Cuchar(3)

const VRNA_DECOMP_PAIR_ML_EXT = Cuchar(23)

const VRNA_DECOMP_PAIR_ML_OUTSIDE = Cuchar(4)

const VRNA_DECOMP_ML_ML_ML = Cuchar(5)

const VRNA_DECOMP_ML_STEM = Cuchar(6)

const VRNA_DECOMP_ML_ML = Cuchar(7)

const VRNA_DECOMP_ML_UP = Cuchar(8)

const VRNA_DECOMP_ML_ML_STEM = Cuchar(9)

const VRNA_DECOMP_ML_COAXIAL = Cuchar(10)

const VRNA_DECOMP_ML_COAXIAL_ENC = Cuchar(11)

const VRNA_DECOMP_EXT_EXT = Cuchar(12)

const VRNA_DECOMP_EXT_UP = Cuchar(13)

const VRNA_DECOMP_EXT_STEM = Cuchar(14)

const VRNA_DECOMP_EXT_EXT_EXT = Cuchar(15)

const VRNA_DECOMP_EXT_STEM_EXT = Cuchar(16)

const VRNA_DECOMP_EXT_STEM_OUTSIDE = Cuchar(17)

const VRNA_DECOMP_EXT_EXT_STEM = Cuchar(18)

const VRNA_DECOMP_EXT_EXT_STEM1 = Cuchar(19)

const VRNA_DECOMP_EXT_STEM_EXT1 = Cuchar(20)

const VRNA_DECOMP_EXT_L = Cuchar(21)

const VRNA_DECOMP_EXT_EXT_L = Cuchar(22)

const VRNA_CONSTRAINT_NO_HEADER = 0

const VRNA_CONSTRAINT_DB = Cuint(16384)

const VRNA_CONSTRAINT_DB_ENFORCE_BP = Cuint(32768)

const VRNA_CONSTRAINT_DB_PIPE = Cuint(65536)

const VRNA_CONSTRAINT_DB_DOT = Cuint(131072)

const VRNA_CONSTRAINT_DB_X = Cuint(262144)

const VRNA_CONSTRAINT_DB_ANG_BRACK = Cuint(524288)

const VRNA_CONSTRAINT_DB_RND_BRACK = Cuint(1048576)

const VRNA_CONSTRAINT_DB_INTRAMOL = Cuint(2097152)

const VRNA_CONSTRAINT_DB_INTERMOL = Cuint(4194304)

const VRNA_CONSTRAINT_DB_GQUAD = Cuint(8388608)

const VRNA_CONSTRAINT_DB_CANONICAL_BP = Cuint(16777216)

const VRNA_CONSTRAINT_DB_WUSS = Cuint(33554432)

const VRNA_CONSTRAINT_DB_DEFAULT = (((((((VRNA_CONSTRAINT_DB | VRNA_CONSTRAINT_DB_PIPE) | VRNA_CONSTRAINT_DB_DOT) | VRNA_CONSTRAINT_DB_X) | VRNA_CONSTRAINT_DB_ANG_BRACK) | VRNA_CONSTRAINT_DB_RND_BRACK) | VRNA_CONSTRAINT_DB_INTRAMOL) | VRNA_CONSTRAINT_DB_INTERMOL) | VRNA_CONSTRAINT_DB_GQUAD

const VRNA_CONSTRAINT_CONTEXT_EXT_LOOP = Cuchar(0x01)

const VRNA_CONSTRAINT_CONTEXT_HP_LOOP = Cuchar(0x02)

const VRNA_CONSTRAINT_CONTEXT_INT_LOOP = Cuchar(0x04)

const VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC = Cuchar(0x08)

const VRNA_CONSTRAINT_CONTEXT_MB_LOOP = Cuchar(0x10)

const VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC = Cuchar(0x20)

const VRNA_CONSTRAINT_CONTEXT_ENFORCE = Cuchar(0x40)

const VRNA_CONSTRAINT_CONTEXT_NO_REMOVE = Cuchar(0x80)

const VRNA_CONSTRAINT_CONTEXT_NONE = Cuchar(0)

const VRNA_CONSTRAINT_CONTEXT_CLOSING_LOOPS = Cuchar(((VRNA_CONSTRAINT_CONTEXT_EXT_LOOP | VRNA_CONSTRAINT_CONTEXT_HP_LOOP) | VRNA_CONSTRAINT_CONTEXT_INT_LOOP) | VRNA_CONSTRAINT_CONTEXT_MB_LOOP)

const VRNA_CONSTRAINT_CONTEXT_ENCLOSED_LOOPS = Cuchar(VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC | VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)

const VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS = Cuchar(VRNA_CONSTRAINT_CONTEXT_CLOSING_LOOPS | VRNA_CONSTRAINT_CONTEXT_ENCLOSED_LOOPS)

const VRNA_CONSTRAINT_WINDOW_UPDATE_5 = Cuint(1)

const VRNA_CONSTRAINT_WINDOW_UPDATE_3 = Cuint(2)

const VRNA_BRACKETS_ALPHA = Cuint(4)

const VRNA_BRACKETS_RND = Cuint(8)

const VRNA_BRACKETS_CLY = Cuint(16)

const VRNA_BRACKETS_ANG = Cuint(32)

const VRNA_BRACKETS_SQR = Cuint(64)

const VRNA_BRACKETS_DEFAULT = ((VRNA_BRACKETS_RND | VRNA_BRACKETS_CLY) | VRNA_BRACKETS_ANG) | VRNA_BRACKETS_SQR

const VRNA_BRACKETS_ANY = (((VRNA_BRACKETS_RND | VRNA_BRACKETS_CLY) | VRNA_BRACKETS_ANG) | VRNA_BRACKETS_SQR) | VRNA_BRACKETS_ALPHA

const VRNA_PLIST_TYPE_BASEPAIR = 0

const VRNA_PLIST_TYPE_GQUAD = 1

const VRNA_PLIST_TYPE_H_MOTIF = 2

const VRNA_PLIST_TYPE_I_MOTIF = 3

const VRNA_PLIST_TYPE_UD_MOTIF = 4

const VRNA_PLIST_TYPE_STACK = 5

const VRNA_STRUCTURE_TREE_HIT = Cuint(1)

const VRNA_STRUCTURE_TREE_SHAPIRO_SHORT = Cuint(2)

const VRNA_STRUCTURE_TREE_SHAPIRO = Cuint(3)

const VRNA_STRUCTURE_TREE_SHAPIRO_EXT = Cuint(4)

const VRNA_STRUCTURE_TREE_SHAPIRO_WEIGHT = Cuint(5)

const VRNA_STRUCTURE_TREE_EXPANDED = Cuint(6)

const VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP = Cuint(1)

const VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP = Cuint(2)

const VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP = Cuint(4)

const VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP = Cuint(8)

const VRNA_UNSTRUCTURED_DOMAIN_MOTIF = Cuint(16)

const VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS = ((VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP | VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP) | VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP) | VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP

const VRNA_MOVESET_INSERTION = 4

const VRNA_MOVESET_DELETION = 8

const VRNA_MOVESET_SHIFT = 16

const VRNA_MOVESET_NO_LP = 32

const VRNA_MOVESET_DEFAULT = VRNA_MOVESET_INSERTION | VRNA_MOVESET_DELETION

const VRNA_VERBOSITY_QUIET = -1

const VRNA_VERBOSITY_DEFAULT = 1

const VRNA_STATUS_MFE_PRE = Cuchar(1)

const VRNA_STATUS_MFE_POST = Cuchar(2)

const VRNA_STATUS_PF_PRE = Cuchar(3)

const VRNA_STATUS_PF_POST = Cuchar(4)

const VRNA_OPTION_DEFAULT = Cuint(0)

const VRNA_OPTION_MFE = Cuint(1)

const VRNA_OPTION_HYBRID = Cuint(4)

const VRNA_OPTION_EVAL_ONLY = Cuint(8)

const VRNA_OPTION_WINDOW = Cuint(16)

const VRNA_NEIGHBOR_CHANGE = 1

const VRNA_NEIGHBOR_INVALID = 2

const VRNA_NEIGHBOR_NEW = 3

const VRNA_PARAMETER_FORMAT_DEFAULT = 0

const VRNA_PBACKTRACK_DEFAULT = 0

const VRNA_PBACKTRACK_NON_REDUNDANT = 1

const VRNA_PLOT_TYPE_SIMPLE = 0

const VRNA_PLOT_TYPE_NAVIEW = 1

const VRNA_PLOT_TYPE_CIRCULAR = 2

const VRNA_PLOT_TYPE_TURTLE = 3

const VRNA_PLOT_TYPE_PUZZLER = 4

const MAXDOS = 1000

const VRNA_UNSORTED = 0

const VRNA_SORT_BY_ENERGY_LEXICOGRAPHIC_ASC = 1

const VRNA_SORT_BY_ENERGY_ASC = 2

# Skipping MacroDefinition: PRIVATE static

const VRNA_INPUT_ERROR = Cuint(1)

const VRNA_INPUT_QUIT = Cuint(2)

const VRNA_INPUT_MISC = Cuint(4)

const VRNA_INPUT_FASTA_HEADER = Cuint(8)

const VRNA_INPUT_SEQUENCE = Cuint(16)

const VRNA_INPUT_CONSTRAINT = Cuint(32)

const VRNA_INPUT_NO_TRUNCATION = Cuint(256)

const VRNA_INPUT_NO_REST = Cuint(512)

const VRNA_INPUT_NO_SPAN = Cuint(1024)

const VRNA_INPUT_NOSKIP_BLANK_LINES = Cuint(2048)

const VRNA_INPUT_BLANK_LINE = Cuint(4096)

const VRNA_INPUT_NOSKIP_COMMENTS = Cuint(128)

const VRNA_INPUT_COMMENT = Cuint(8192)

# exports
const PREFIXES = ["VRNA_", "vrna_"]
for name in names(@__MODULE__; all=true), prefix in PREFIXES
    if startswith(string(name), prefix)
        @eval export $name
    end
end

end # module
