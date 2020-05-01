#include "declarations.h"

double Recount_Lip(const int comp_num, const int sign, const double * lipids, const double ultimate)
{
    if(comp_num == 1)
        return ultimate * 100 / lipids[0];
    //----------------------------------
    if(comp_num == 2 && sign == 1)
        return ultimate * 100 / lipids[0];
    if(comp_num == 2 && sign == 2)
        return ultimate * 100 / lipids[1];
    //----------------------------------
    if(comp_num == 3 && sign == 1)
        return ultimate * 100 / lipids[0];
    if(comp_num == 3 && sign == 2)
        return ultimate * 100 / lipids[1];
    if(comp_num == 3 && sign == 3)
        return ultimate * 100 / lipids[2];
    //----------------------------------
    if(comp_num == 4 && sign == 1)
        return ultimate * 100 / lipids[0];
    if(comp_num == 4 && sign == 2)
        return ultimate * 100 / lipids[1];
    if(comp_num == 4 && sign == 3)
        return ultimate * 100 / lipids[2];
    if(comp_num == 4 && sign == 4)
        return ultimate * 100 / lipids[3];
    //----------------------------------
    if(comp_num == 5 && sign == 1)
        return ultimate * 100 / lipids[0];
    if(comp_num == 5 && sign == 2)
        return ultimate * 100 / lipids[1];
    if(comp_num == 5 && sign == 3)
        return ultimate * 100 / lipids[2];
    if(comp_num == 5 && sign == 4)
        return ultimate * 100 / lipids[3];
    if(comp_num == 5 && sign == 5)
        return ultimate * 100 / lipids[4];
}
//======================================================================================
double Ratio_Calc(const int comp_num, const double* prop, const double recount_lip[][LIP_PROP], const int i)
{
    if(comp_num == 1)
        return recount_lip[0][i];
    if(comp_num == 2)
        return (prop[0] * recount_lip[0][i]) + (prop[1] * recount_lip[1][i]);
    if(comp_num == 3)
        return (prop[0] * recount_lip[0][i]) + (prop[1] * recount_lip[1][i]) + (prop[2] * recount_lip[2][i]);
    if(comp_num == 4)
        return (prop[0] * recount_lip[0][i]) + (prop[1] * recount_lip[1][i]) + (prop[2] * recount_lip[2][i])
               + (prop[3] * recount_lip[3][i]);
    if(comp_num == 5)
        return (prop[0] * recount_lip[0][i]) + (prop[1] * recount_lip[1][i]) + (prop[2] * recount_lip[2][i])
               + (prop[3] * recount_lip[3][i]) + (prop[4] * recount_lip[4][i]);
}
//======================================================================================
double Lip_Balance_Ratio(const double fao_voz2008, const double recount_lip)
{
    if(recount_lip <= fao_voz2008)
        return recount_lip / fao_voz2008;
    else
        return fao_voz2008 / recount_lip;
}
//======================================================================================
void Fattyacid_Compliance(const double * lip_balance_ratio, double* result1, double* result2)
{
    double multi1, multi2;
    multi1 = multi2 = 1.0;
    for(int i = 0; i < 3; i++)
        multi1 *= lip_balance_ratio[i];
    *result1 = pow(multi1, 1.0 / 3.0);
    for(int i = 0; i < LIP_PROP; i++)
        multi2 *= lip_balance_ratio[i];
    *result2 = pow(multi2, 1.0 / 5.0);
}
