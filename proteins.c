#include "declarations.h"

/**
 * @file     proteins.c
 * @Author   Daniil Kuprianov (vaedermakar@mail.ru)
 * @date     September, 2019
 * @brief    proteins definitions of the program
 */

// переменная sign даёт понять функции, какой именно компонент требуется
// пересчитать во время её текущего вызова
double Recount(int comp_num, int sign, const double * protein, const double comp)
{
    if(comp_num == 1)
        return comp * 100 / protein[0];
    //----------------------------------
    if(comp_num == 2 && sign == 1)
        return comp * 100 / protein[0];
    if(comp_num == 2 && sign == 2)
        return comp * 100 / protein[1];
    //----------------------------------
    if(comp_num == 3 && sign == 1)
        return comp * 100 / protein[0];
    if(comp_num == 3 && sign == 2)
        return comp * 100 / protein[1];
    if(comp_num == 3 && sign == 3)
        return comp * 100 / protein[2];
    //----------------------------------
    if(comp_num == 4 && sign == 1)
        return comp * 100 / protein[0];
    if(comp_num == 4 && sign == 2)
        return comp * 100 / protein[1];
    if(comp_num == 4 && sign == 3)
        return comp * 100 / protein[2];
    if(comp_num == 4 && sign == 4)
        return comp * 100 / protein[3];
    //----------------------------------
    if(comp_num == 5 && sign == 1)
        return comp * 100 / protein[0];
    if(comp_num == 5 && sign == 2)
        return comp * 100 / protein[1];
    if(comp_num == 5 && sign == 3)
        return comp * 100 / protein[2];
    if(comp_num == 5 && sign == 4)
        return comp * 100 / protein[3];
    if(comp_num == 5 && sign == 5)
        return comp * 100 / protein[4];
}
//======================================================================================
double AKP(int comp_num, const double * prop, const double recount[][AMI], int i)
{
    if(comp_num == 1)
        return recount[0][i];
    if(comp_num == 2)
        return (prop[0] * recount[0][i]) + (prop[1] * recount[1][i]);
    if(comp_num == 3)
        return (prop[0] * recount[0][i]) + (prop[1] * recount[1][i]) + (prop[2] * recount[2][i]);
    if(comp_num == 4)
        return (prop[0] * recount[0][i]) + (prop[1] * recount[1][i]) + (prop[2] * recount[2][i])
               + (prop[3] * recount[3][i]);
    if(comp_num == 5)
        return (prop[0] * recount[0][i]) + (prop[1] * recount[1][i]) + (prop[2] * recount[2][i])
               + (prop[3] * recount[3][i]) + (prop[4] * recount[4][i]);
}
//======================================================================================
double Aminoacidskor(const double akp, const double fao_voz2007)
{
    return akp / fao_voz2007 * 100;
}
//======================================================================================
double Koef_Ration(double min_c, const double aminoacidskor)
{
    return min_c / aminoacidskor;
}
//======================================================================================
double Aminoacidskor_Sum(const double * aminoacidskor)
{
    double sum = 0.0;
    for(int i = 0; i < AMI; i++)
    {
        sum += *aminoacidskor;
        aminoacidskor++;
    }
    return sum;
}
//======================================================================================
double KRAS(const double aminoacidskor_sum, double min_c)
{
    return (aminoacidskor_sum - (min_c * AMI)) / AMI;
}
//======================================================================================
double Biological_Value(double kras)
{
    return 100 - kras;
}
//======================================================================================
double Amino_Acid_Comp_Ratio_Coef(const double * koef_ration, const double * akp)
{
    double numerator = 0.0; // числитель в формуле
    double denominator = 0.0; // знаменатель в формуле
    for(int i = 0; i < AMI; i++)
    {
        numerator += koef_ration[i] * akp[i];
        denominator += akp[i];
    }
    return numerator / denominator;
}
//======================================================================================
double Comparable_Redundancy_Ratio(const double * akp, double min_c, const double * fao_voz2007)
{
    double numerator = 0.0;
    for(int i = 0; i < AMI; i++)
    {
        numerator += akp[i] - min_c / 100 * fao_voz2007[i];
    }
    return numerator / (min_c / 100);
}
//======================================================================================
double Fatty_Acid_Per_100g(const double fao_voz2007, const double akp)
{
    if(fao_voz2007 <= akp)
        return fao_voz2007 / akp;
    else
        return akp / fao_voz2007;
}
//======================================================================================
double Balance_Index(const double * fatty_acid_per_100g)
{
    double temp = 1.0;
    for(int i = 0; i < AMI; i++)
    {
        temp *= fatty_acid_per_100g[i];
    }
    return pow(temp, 1.0/9.0);
}
//======================================================================================
double Balance_Index_General(const double balance_index, const double biological_value, const double amino_acid_comp_ratio_coef)
{
    double temp = (balance_index / 1) * (biological_value / 100) * (amino_acid_comp_ratio_coef / 1);
    return pow(temp, 1.0/3.0);
}
