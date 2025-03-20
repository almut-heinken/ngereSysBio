function model=defineMediumDMEM(model)
% define constraints for DMEM (Dulbecco's Modified Eagle Medium) with high glucose

% https://www.sigmaaldrich.com/FR/fr/product/sigma/d5796

% based on previous DMEM definition from PMID:31042485

model=changeRxnBounds(model,model.rxns(find(strncmp(model.rxns,'EX_',3))),0,'l');

%% amino acids
model=changeRxnBounds(model,'EX_gly[e]',-1,'l');
model=changeRxnBounds(model,'EX_arg_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_cys_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_gln_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_his_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_ile_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_leu_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_lys_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_met_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_phe_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_pro_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_ser_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_thr_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_trp_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_tyr_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_val_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_ala_L[e]',-1,'l');
model=changeRxnBounds(model,'EX_glc_D[e]',-5,'l');
%% other
model=changeRxnBounds(model,'EX_chol[e]',-1,'l');
model=changeRxnBounds(model,'EX_pyr[e]',-1,'l');
model=changeRxnBounds(model,'EX_lnlc[e]',-1,'l');
model=changeRxnBounds(model,'EX_lipoate[e]',-1,'l');

%% ions and vitamins: no exact composition given in the experimental medium
% I assume the ones required by the reconstructions
% ions
model=changeRxnBounds(model,'EX_ca2[e]',-1,'l');
model=changeRxnBounds(model,'EX_cl[e]',-1,'l');
model=changeRxnBounds(model,'EX_so4[e]',-1,'l');
model=changeRxnBounds(model,'EX_h2s[e]',-1,'l');
model=changeRxnBounds(model,'EX_cobalt2[e]',-1,'l');
model=changeRxnBounds(model,'EX_cu2[e]',-1,'l');
model=changeRxnBounds(model,'EX_fe2[e]',-1,'l');
model=changeRxnBounds(model,'EX_fe3[e]',-1,'l');
model=changeRxnBounds(model,'EX_k[e]',-1,'l');
model=changeRxnBounds(model,'EX_mg2[e]',-1,'l');
model=changeRxnBounds(model,'EX_mn2[e]',-1,'l');
model=changeRxnBounds(model,'EX_zn2[e]',-1,'l');
model=changeRxnBounds(model,'EX_pi[e]',-10,'l');
model=changeRxnBounds(model,'EX_h2o[e]',-100,'l');
model=changeRxnBounds(model,'EX_co2[e]',-100,'l');
model=changeRxnBounds(model,'EX_hco3[e]',-100,'l');

%% vitamins
model=changeRxnBounds(model,'EX_aqcobal[e]',-1,'l');
model=changeRxnBounds(model,'EX_ccbl[e]',-1,'l');
model=changeRxnBounds(model,'EX_hoxocbl[e]',-1,'l');
model=changeRxnBounds(model,'EX_C06453[e]',-1,'l');
model=changeRxnBounds(model,'EX_btn[e]',-1,'l');
model=changeRxnBounds(model,'EX_fol[e]',-1,'l');
model=changeRxnBounds(model,'EX_thf[e]',-1,'l');
model=changeRxnBounds(model,'EX_ncam[e]',-1,'l');
model=changeRxnBounds(model,'EX_pnto_R[e]',-1,'l');
model=changeRxnBounds(model,'EX_pydxn[e]',-1,'l');
model=changeRxnBounds(model,'EX_ribflv[e]',-1,'l');
model=changeRxnBounds(model,'EX_thm[e]',-1,'l');
model=changeRxnBounds(model,'EX_inost[e]',-1,'l');

% oxygen
model=changeRxnBounds(model,'EX_o2[e]',-100,'l');

end