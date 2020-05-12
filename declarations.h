#pragma once

/**
 * @file     declarations.h
 * @Author   Daniil Kuprianov (vaedermakar@mail.ru)
 * @date     September, 2019
 * @brief    declarations of the program
 */

#define AMI 9 // количество аминокислот
#define LIP_PROP 5 // липидные хар-ки (аналог аминокислот)
#define MAX_COMP 5 // максимальное количество компонентов

//======================================================================================
/*БЕЛКИ*/
// пересчёт аминокислот со 100г компонента на 100г белка
double Recount(const int comp_num, const int sign, const double * protein, const double comp);
// аминокислота в продукте
double AKP(const int comp_num, const double * prop, const double recount[][AMI], const int i);
// аминокислотный скор
double Aminoacidskor(const double akp, const double fao_voz2007);
// коэффициент рациональности
double Koef_Ration(const double min_c, const double aminoacidskor);
// сумма аминокислотных скоров для всех аминокислот
double Aminoacidskor_Sum(const double * aminoacidskor);
// коэффициент разбалансированности
double KRAS(const double aminoacidskor_sum, const double min_c);
// биологическая ценность
double Biological_Value(const double kras);
// коэффициент рациональности аминокислотного состава
double Amino_Acid_Comp_Ratio_Coef(const double * koef_ration, const double * akp);
// показатель сопоставимой избыточности
double Comparable_Redundancy_Ratio(const double * akp, const double min_c, const double * fao_voz2007);
// продукт с учётом пропорций (жирная кислота на 100г продукта, аналог АКП)
double Fatty_Acid_Per_100g(const double fao_voz2007, const double akp);
// коэффициент сбалансированности
double Balance_Index(const double * fatty_acid_per_100g);
// коэффициент сбалансированности всей белковой составляющей
double Balance_Index_General(const double balance_index, const double biological_value, const double amino_acid_comp_ratio_coef);
//======================================================================================
/*ЛИПИДЫ*/
// пересчёт на 100г жира всех липидных входных данных
double Recount_Lip(const int comp_num, const int sign, const double * lipids, const double ultimate);
// коэффициент рациональности жирной кислоты
double Ratio_Calc(const int comp_num, const double * prop, const double recount_lip[][LIP_PROP], const int i);
// коэффициент сбалансированности
double Lip_Balance_Ratio(const double fao_voz2008, const double recount_lip);
// коэффициент жирнокислотного соответствия
void Fattyacid_Compliance(const double * lip_balance_ratio, double* result1, double* result2);
