/*Программатор Куприянов Д.К., осень 2019 года.
 Программа производит оценку белковой  и липидной составляющих
 продуктов питания.*/
 /*
 1.Вводим кол-во компонентов
 2.Вводим кол-ва белка на 100г каждого компонента (после ввода кол-ва компонентов) (n в формуле)(от 1 до 5, значение от 0 до 100)
 3.Вводим пропорции компонентов
 4.Вводим аминокислоты 9 раз для каждого компонента (Xётая в формуле)
 5.Пересчитываем аминокислоты с продукта на белок
 4.Вычисляем АКП, С, Аётую, крас, БЦ, R, G
*/
#include <stdio.h>
#include <windows.h>
#include <math.h>
#include <stdlib.h>

#define AMI 9 // количество аминокислот
#define LIP_PROP 5 // липидные хар-ки (аналог аминокислот)
#define MAX_COMP 5 // максимальное количество компонентов

// данная версия программы будет запускаться только на флешке
// и предназначена для тестовой демонстрации
// #define FLASH_TEST
//======================================================================================
/*БЕЛКИ*/
// пересчёт аминокислот со 100г компонента на 100г белка
double Recount(int comp_num, int sign, double * protein, double comp);
// аминокислота в продукте
double AKP(int comp_num, double * prop, double recount[][AMI], int i);
// аминокислотный скор
double Aminoacidskor(double akp, const double fao_voz2007);
// коэффициент рациональности
double Koef_Ration(double min_c, double c);
// сумма аминокислотных скоров для всех аминокислот
double C_SUM(double * c);
// коэффициент разбалансированности
double KRAS(double c_sum, double min_c);
// биологическая ценность
double BC(double kras);
// коэффициент рациональности аминокислотного состава
double R(double * aj, double * akp);
// показатель сопоставимой избыточности
double G(double * akp, double min_c, const double * fao_voz2007);
// продукт с учётом пропорций (жирная кислота на 100г продукта, аналог АКП)
double K(const double fao_voz2007, double akp);
// коэффициент сбалансированности
double S(double * k);
// коэффициент сбалансированности всей белковой составляющей
double K_general(double s, double bc, double r);
//======================================================================================
/*ЛИПИДЫ*/
// пересчёт на 100г жира всех липидных входных данных
double Recount_lip(int comp_num, int sign, double * lipids, double ultimate);
// коэффициент рациональности жирной кислоты
double Ratio_calc(int comp_num, double * prop, double recount_lip[][LIP_PROP], int i);
// коэффициент сбалансированности
double Lip_balance_ratio(double fao_voz2008, double recount_lip);
// коэффициент жирнокислотного соответствия
void Fattyacid_compliance(double * Lip_balance_ratio, double* result1, double* result2);

int main(void)
{
	SetConsoleCP(1251);SetConsoleOutputCP(1251);
	#ifdef FLASH_TEST
	int is_flash = GetDriveType(NULL); // NULL подразумевает что проверяется текущий носитель
	if(is_flash != 2)
    {
        printf("Ошибка доступа!");
        exit(EXIT_FAILURE);
    }
    #endif // FLASH_TEST
    int ans = 0;
    system("color 70");
    puts("По всем вопросам - vaedermakar@mail.ru.\n");
    puts("*Программа использует ФАО ВОЗ 2007 для вычислений белковой составляющей\n"
        "ФАО ВОЗ 2008 для вычислений липидной составляющей.");
    do{
    /*ДЛЯ ХРАНЕНИЯ ДАННЫХ*/
    /*ОБЩЕЕ*/
    int comp_num = 0; // кол-во компонентов
    /*БЕЛКИ*/
    const double fao_voz2007[AMI] = {3.9, 1.5, 3.0, 5.9, 4.5, 2.2, 2.3, 0.6, 3.8};
    double protein[MAX_COMP] = {0.0}; // кол-во белка на 100г каждого компонента
	double comp[MAX_COMP][AMI] = {0.0}; // для аминокислот
	/*ЛИПИДЫ*/
    const double fao_voz2008[MAX_COMP] = {33.33, 33.33, 33.33, 6.67, 26.67};
    double lipids[MAX_COMP] = {0.0}; // кол-во липидов на 100г каждого компонента
    /*Здесь один массив отвечает за ненасыщенные жирные кислоты,
    мононенасыщенные жирные кислоты, полиненасыщенные жирные кислоты, омега3 и омега6.
    Это сделано для того чтобы передать один двумерный массив в функцию вместо 5 обычных.
    [0] - насыщенные жирные кислоты
    [1] - мононенасыщенные жирные кислоты
    [2] - полиненасыщенные жирные кислоты
    [3] - омега3
    [4] - омега6*/
    double ultimate[MAX_COMP][LIP_PROP] = {0.0};

	/*ДЛЯ РЕЗУЛЬТАТОВ ВЫЧИСЛЕНИЙ*/
	/*БЕЛКИ*/
	double akp[AMI] = {0.0};
	double c[AMI] = {0.0};
	double aj[AMI] = {0.0};
	double k[AMI] = {0.0};
	double prop[MAX_COMP] = {0.0}; // для пропорций, если компонентов > 1
	double recount[MAX_COMP][AMI] = {0.0}; // для пересчёта аминокилот
    /*коэффициент разбалансированности, биологическая ценность,
    коэффициент рациональности аминокислотного состава,
    индекс сбалансированности, показатель сопоставимой избыточности и
    коэффициент сбалансированности всей белковой составляющей(вычисляются 1 раз в конце)*/
    double kras = 0.0;
    double bc = 0.0;
    double r = 0.0;
    double g = 0.0;
    double s = 0.0;
    double k_general = 0.0;
    int i; // счётчик, номер строки с аминокислотой

    /*ЛИПИДЫ*/
    double recount_lip[MAX_COMP][LIP_PROP] = {0.0};
    double ratio_calc[LIP_PROP] = {0.0};
	double lip_balance_ratio[LIP_PROP] = {0.0};
	double result1 = 0.0;
	double result2 = 0.0;
    //======================================================================================
    do{
    puts("Введите количество компонентов продукта (до 5):");
    while(scanf(" %d", &comp_num) != 1)
    {
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    putchar('\n');
    // если компонентов 0 или меньше либо больше 5
    if(comp_num <= 0 || comp_num > MAX_COMP)
        printf("Компонентов не может быть %d, попробуйте ещё раз.\n", comp_num);
    }while(comp_num <= 0 || comp_num > MAX_COMP);
    //======================================================================================
    do{
    if(comp_num == 1)
    {
        prop[0] = 1.0;
    }
    if(comp_num > 1)
    {
        puts("Введите пропорции компонентов(таким образом - 0.20 0.60 0.20):");
        if(comp_num == 2)
        {
            while(scanf(" %lf %lf", &prop[0], &prop[1]) != 2)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
        }
        if(comp_num == 3)
        {
            while(scanf(" %lf %lf %lf", &prop[0], &prop[1], &prop[2]) != 3)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
        }
        if(comp_num == 4)
        {
            while(scanf(" %lf %lf %lf %lf", &prop[0], &prop[1], &prop[2], &prop[3]) != 4)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
        }
        if(comp_num == 5)
        {
            while(scanf(" %lf %lf %lf %lf %lf", &prop[0], &prop[1], &prop[2], &prop[3], &prop[4]) != 5)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
        }
        putchar('\n');
    }
    // модуль разности сравниваем с дельтой
    // если компонентов больше 1 и не больше 5, сумма элементов должна равняться 1.0(100%)
    if(fabs((prop[0]+prop[1]+prop[2]+prop[3]+prop[4])-(1.0)) >= 0.000000000001)
        puts("Ошибка ввода. Сумма пропорций компонентов должна равняться 1.0.");
	}while(fabs((prop[0]+prop[1]+prop[2]+prop[3]+prop[4])-(1.0)) >= 0.000000000001);
    //======================================================================================
    if(comp_num > 1)
        puts("Введите кол-во белка на 100г каждого компонента:");
    if(comp_num == 1)
    {
        puts("Введите количество белка на 100г компонента:");
        while(scanf(" %lf", &protein[0]) != 1)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    if(comp_num == 2)
    {
        while(scanf(" %lf %lf", &protein[0], &protein[1]) != 2
            || protein[0] < 0.0 || protein[0] > 100.0 || protein[1] < 0.0 || protein[1] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    if(comp_num == 3)
    {
        while(scanf(" %lf %lf %lf", &protein[0], &protein[1], &protein[2]) != 3
            || protein[0] < 0.0 || protein[0] > 100.0 || protein[1] < 0.0 || protein[1] > 100.0
            || protein[2] < 0.0 || protein[2] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    if(comp_num == 4)
    {
        while(scanf(" %lf %lf %lf %lf", &protein[0], &protein[1], &protein[2], &protein[3]) != 4
            || protein[0] < 0.0 || protein[0] > 100.0 || protein[1] < 0.0 || protein[1] > 100.0
            || protein[2] < 0.0 || protein[2] > 100.0 || protein[3] < 0.0 || protein[3] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    if(comp_num == 5)
    {
        while(scanf(" %lf %lf %lf %lf %lf", &protein[0], &protein[1], &protein[2], &protein[3], &protein[4]) != 5
            || protein[0] < 0.0 || protein[0] > 100.0 || protein[1] < 0.0 || protein[1] > 100.0
            || protein[2] < 0.0 || protein[2] > 100.0 || protein[3] < 0.0 || protein[3] > 100.0
            || protein[4] < 0.0 || protein[4] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    putchar('\n');
    //======================================================================================
    if(comp_num > 1)
        puts("Введите аминокислоты в граммах для каждого компонента.");
    if(comp_num == 1)
        puts("Введите аминокислоты в граммах.");
    for(i = 0; i < AMI; i++)
    {
        if(i == 0 && comp_num > 1)
            puts("Введите Валин для каждого компонента:");
        if(i == 1 && comp_num > 1)
            puts("Введите Гистидин для каждого компонента:");
        if(i == 2 && comp_num > 1)
           puts("Введите Изолейцин для каждого компонента:");
        if(i == 3 && comp_num > 1)
           puts("Введите Лейцин для каждого компонента:");
        if(i == 4 && comp_num > 1)
            puts("Введите Лизин для каждого компонента:");
        if(i == 5 && comp_num > 1)
            puts("Введите Метионин+Цистеин для каждого компонента:");
        if(i == 6 && comp_num > 1)
            puts("Введите Треонин для каждого компонента:");
        if(i == 7 && comp_num > 1)
            puts("Введите Триптофан для каждого компонента:");
        if(i == 8 && comp_num > 1)
           puts("Введите Фенилаланин+Тирозин для каждого компонента:");

        if(i == 0 && comp_num == 1)
            puts("Введите Валин:");
        if(i == 1 && comp_num == 1)
            puts("Введите Гистидин:");
        if(i == 2 && comp_num == 1)
            puts("Введите Изолейцин:");
        if(i == 3 && comp_num == 1)
            puts("Введите Лейцин:");
        if(i == 4 && comp_num == 1)
            puts("Введите Лизин:");
        if(i == 5 && comp_num == 1)
            puts("Введите Метионин+Цистеин:");
        if(i == 6 && comp_num == 1)
            puts("Введите Треонин:");
        if(i == 7 && comp_num == 1)
            puts("Введите Триптофан:");
        if(i == 8 && comp_num == 1)
            puts("Введите Фенилаланин+Тирозин:");

        if(comp_num == 1)
        {
            while(scanf(" %lf", &comp[0][i]) != 1 || comp[0][i] <= 0.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 2)
        {
            while(scanf(" %lf %lf", &comp[0][i], &comp[1][i]) != 2
                  || comp[0][i] <= 0.0 || comp[1][i] <= 0.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 3)
        {
            while(scanf(" %lf %lf %lf", &comp[0][i], &comp[1][i], &comp[2][i]) != 3
                  || comp[0][i] <= 0.0 || comp[1][i] <= 0.0 || comp[2][i] <= 0.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
    	if(comp_num == 4)
        {
            while(scanf(" %lf %lf %lf %lf", &comp[0][i], &comp[1][i], &comp[2][i], &comp[3][i]) != 4
                  || comp[0][i] <= 0.0 || comp[1][i] <= 0.0 || comp[2][i] <= 0.0 || comp[3][i] <= 0.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 5)
        {
            while(scanf(" %lf %lf %lf %lf %lf", &comp[0][i], &comp[1][i], &comp[2][i], &comp[3][i], &comp[4][i]) != 5
                  || comp[0][i] <= 0.0 || comp[1][i] <= 0.0 || comp[2][i] <= 0.0 || comp[3][i] <= 0.0 || comp[4][i] <= 0.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
    }
    //======================================================================================
    // ввод жиров на 100г каждого компонента для пересчета
    if(comp_num > 1)
        puts("Введите кол-во липидов на 100г каждого компонента:");
    if(comp_num == 1)
    {
        puts("Введите кол-во липидов на 100г компонента:");
        while(scanf(" %lf", &lipids[0]) != 1 || lipids[0] < 0.0 || lipids[0] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    if(comp_num == 2)
    {
        while(scanf(" %lf %lf", &lipids[0], &lipids[1]) != 2
            || lipids[0] < 0.0 || lipids[0] > 100.0 || lipids[1] < 0.0 || lipids[1] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    if(comp_num == 3)
    {
        while(scanf(" %lf %lf %lf", &lipids[0], &lipids[1], &lipids[2]) != 3
            || lipids[0] < 0.0 || lipids[0] > 100.0 || lipids[1] < 0.0 || lipids[1] > 100.0
            || lipids[2] < 0.0 || lipids[2] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    if(comp_num == 4)
    {
        while(scanf(" %lf %lf %lf %lf", &lipids[0], &lipids[1], &lipids[2], &lipids[3]) != 4
            || lipids[0] < 0.0 || lipids[0] > 100.0 || lipids[1] < 0.0 || lipids[1] > 100.0
            || lipids[2] < 0.0 || lipids[2] > 100.0 || lipids[3] < 0.0 || lipids[3] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    if(comp_num == 5)
    {
        while(scanf(" %lf %lf %lf %lf %lf", &lipids[0], &lipids[1], &lipids[2], &lipids[3], &lipids[4]) != 5
            || lipids[0] < 0.0 || lipids[0] > 100.0 || lipids[1] < 0.0 || lipids[1] > 100.0
            || lipids[2] < 0.0 || lipids[2] > 100.0 || lipids[3] < 0.0 || lipids[3] > 100.0
            || lipids[4] < 0.0 || lipids[4] > 100.0)
        {
            puts("Некорректный ввод.");
            fflush(stdin);
        }
    }
    putchar('\n');
    //======================================================================================
   /*ввод:
    [0] - насыщенные жирные кислоты
    [1] - мононенасыщенные жирные кислоты
    [2] - полиненасыщенные жирные кислоты
    [3] - омега3
    [4] - омега6*/
    for(i = 0; i < LIP_PROP; i++)
    {
        if(comp_num > 1 && i == 0)
            puts("Введите кол-во насыщенных жирных кислот на 100г каждого компонента:");
        if(comp_num > 1 && i == 1)
            puts("Введите кол-во мононенасыщенных жирных кислот на 100г каждого компонента:");
        if(comp_num > 1 && i == 2)
            puts("Введите кол-во полиненасыщенных жирных кислот на 100г каждого компонента:");
        if(comp_num > 1 && i == 3)
            puts("Введите омега-3 на 100г каждого компонента:");
        if(comp_num > 1 && i == 4)
            puts("Введите омега-6 на 100г каждого компонента:");
        if(comp_num == 1)
        {
            if(i == 0)
                puts("Введите кол-во насыщенных жирных кислот на 100г компонента:");
            if(i == 1)
                puts("Введите кол-во мононенасыщенных жирных кислот на 100г компонента:");
            if(i == 2)
                puts("Введите кол-во полиненасыщенных жирных кислот на 100г компонента:");
            if(i == 3)
                puts("Введите кол-во омега-3 на 100г компонента:");
            if(i == 4)
                puts("Введите кол-во омега-6 на 100г компонента:");

            while(scanf(" %lf", &ultimate[0][i]) != 1 || ultimate[0][i] < 0.0 || ultimate[0][i] > 100.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 2)
        {
            while(scanf(" %lf %lf", &ultimate[0][i], &ultimate[1][i]) != 2
                || ultimate[0][i] < 0.0 || ultimate[0][i] > 100.0 || ultimate[1][i] < 0.0 || ultimate[1][i] > 100.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 3)
        {
            while(scanf(" %lf %lf %lf", &ultimate[0][i], &ultimate[1][i], &ultimate[2][i]) != 3
                || ultimate[0][i] < 0.0 || ultimate[0][i] > 100.0 || ultimate[1][i] < 0.0 || ultimate[1][i] > 100.0
                || ultimate[2][i] < 0.0 || ultimate[2][i] > 100.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 4)
        {
            while(scanf(" %lf %lf %lf %lf", &ultimate[0][i], &ultimate[1][i], &ultimate[2][i], &ultimate[3][i]) != 4
                || ultimate[0][i] < 0.0 || ultimate[0][i] > 100.0 || ultimate[1][i] < 0.0 || ultimate[1][i] > 100.0
                || ultimate[2][i] < 0.0 || ultimate[2][i] > 100.0 || ultimate[3][i] < 0.0 || ultimate[3][i] > 100.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 5)
        {
            while(scanf(" %lf %lf %lf %lf %lf", &ultimate[0][i], &ultimate[1][i], &ultimate[2][i], &ultimate[3][i], &ultimate[4][i]) != 5
                || ultimate[0][i] < 0.0 || ultimate[0][i] > 100.0 || ultimate[1][i] < 0.0 || ultimate[1][i] > 100.0
                || ultimate[2][i] < 0.0 || ultimate[2][i] > 100.0 || ultimate[3][i] < 0.0 || ultimate[3][i] > 100.0
                || ultimate[4][i] < 0.0 || ultimate[4][i] > 100.0)
            {
                puts("Некорректный ввод.");
                fflush(stdin);
            }
            putchar('\n');
        }
    }
    //======================================================================================
    /*ВЫЧИСЛЕНИЯ БЕЛКОВ*/
    for(i = 0; i < AMI; i++)
    {
        if(comp_num == 1)
            recount[0][i] = Recount(comp_num, 1, protein, comp[0][i]);
        if(comp_num == 2)
        {
            recount[0][i] = Recount(comp_num, 1, protein, comp[0][i]);
            recount[1][i] = Recount(comp_num, 2, protein, comp[1][i]);
        }
        if(comp_num == 3)
        {
            recount[0][i] = Recount(comp_num, 1, protein, comp[0][i]);
            recount[1][i] = Recount(comp_num, 2, protein, comp[1][i]);
            recount[2][i] = Recount(comp_num, 3, protein, comp[2][i]);
        }
        if(comp_num == 4)
        {
            recount[0][i] = Recount(comp_num, 1, protein, comp[0][i]);
            recount[1][i] = Recount(comp_num, 2, protein, comp[1][i]);
            recount[2][i] = Recount(comp_num, 3, protein, comp[2][i]);
            recount[3][i] = Recount(comp_num, 4, protein, comp[3][i]);
        }
        if(comp_num == 5)
        {
            recount[0][i] = Recount(comp_num, 1, protein, comp[0][i]);
            recount[1][i] = Recount(comp_num, 2, protein, comp[1][i]);
            recount[2][i] = Recount(comp_num, 3, protein, comp[2][i]);
            recount[3][i] = Recount(comp_num, 4, protein, comp[3][i]);
            recount[4][i] = Recount(comp_num, 5, protein, comp[4][i]);
        }

		akp[i] = AKP(comp_num, prop, recount, i);
    	c[i] = Aminoacidskor(akp[i], fao_voz2007[i]);
        k[i] = K(fao_voz2007[i], akp[i]);
    }

    // в конце находим минимальный 'с' из всех, считаем aj для каждой аминокислоты,
    // коэффициент рациональности аминокислотного состава и показатель сопоставимой избыточности
    double min_c = 999.0;
    for(i = 0; i < AMI; i++)
    {
    	if(c[i] < min_c)
    		min_c = c[i];
    }
    for(i = 0; i < AMI; i++)
    {
    	aj[i] = Koef_Ration(min_c, c[i]);
    }
    //======================================================================================
    /*ВЫЧИСЛЕНИЯ ЛИПИДОВ*/
    // здесь i это номер "свойства" липидов
    for(i = 0; i < LIP_PROP; i++)
    {
        if(comp_num == 1)
            recount_lip[0][i] = Recount_lip(comp_num, 1, lipids, ultimate[0][i]);
        if(comp_num == 2)
        {
            recount_lip[0][i] = Recount_lip(comp_num, 1, lipids, ultimate[0][i]);
            recount_lip[1][i] = Recount_lip(comp_num, 2, lipids, ultimate[1][i]);
        }
        if(comp_num == 3)
        {
            recount_lip[0][i] = Recount_lip(comp_num, 1, lipids, ultimate[0][i]);
            recount_lip[1][i] = Recount_lip(comp_num, 2, lipids, ultimate[1][i]);
            recount_lip[2][i] = Recount_lip(comp_num, 3, lipids, ultimate[2][i]);
        }
        if(comp_num == 4)
        {
            recount_lip[0][i] = Recount_lip(comp_num, 1, lipids, ultimate[0][i]);
            recount_lip[1][i] = Recount_lip(comp_num, 2, lipids, ultimate[1][i]);
            recount_lip[2][i] = Recount_lip(comp_num, 3, lipids, ultimate[2][i]);
            recount_lip[3][i] = Recount_lip(comp_num, 4, lipids, ultimate[3][i]);
        }
        if(comp_num == 5)
        {
            recount_lip[0][i] = Recount_lip(comp_num, 1, lipids, ultimate[0][i]);
            recount_lip[1][i] = Recount_lip(comp_num, 2, lipids, ultimate[1][i]);
            recount_lip[2][i] = Recount_lip(comp_num, 3, lipids, ultimate[2][i]);
            recount_lip[3][i] = Recount_lip(comp_num, 4, lipids, ultimate[3][i]);
            recount_lip[4][i] = Recount_lip(comp_num, 5, lipids, ultimate[4][i]);
        }
        ratio_calc[i] = Ratio_calc(comp_num, prop, recount_lip, i);
        lip_balance_ratio[i] = Lip_balance_ratio(fao_voz2008[i], ratio_calc[i]);
        Fattyacid_compliance(lip_balance_ratio, &result1, &result2);
    }


    //======================================================================================
    /*ВЫВОД АМИНОКИСЛОТ*/
    puts("=================================================");
    for(i = 0; i < AMI; i++)
    {
        if(i == 0)
            puts("Валин:\n");
        if(i == 1)
            puts("Гистидин:\n");
        if(i == 2)
            puts("Изолейцин:\n");
        if(i == 3)
            puts("Лейцин:\n");
        if(i == 4)
            puts("Лизин:\n");
        if(i == 5)
            puts("Метионин+Цистеин:\n");
        if(i == 6)
            puts("Треонин:\n");
        if(i == 7)
            puts("Триптофан:\n");
        if(i == 8)
            puts("Фенилаланин+Тирозин:\n");

    	if(comp_num == 1)
            printf("Аминокислота на 100г белка = %.2f\n", recount[0][i]);
        if(comp_num == 2)
        {
            printf("Аминокислота на 100г белка в 1 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 2 компоненте = %.2f\n",
                    recount[0][i], recount[1][i]);
        }
        if(comp_num == 3)
        {
            printf("Аминокислота на 100г белка в 1 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 2 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 3 компоненте = %.2f\n",
                   recount[0][i], recount[1][i], recount[2][i]);
        }
        if(comp_num == 4)
        {
            printf("Аминокислота на 100г белка в 1 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 2 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 3 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 4 компоненте = %.2f\n",
                   recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
        }
        if(comp_num == 5)
        {
            printf("Аминокислота на 100г белка в 1 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 2 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 3 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 4 компоненте = %.2f\n"
                   "Аминокислота на 100г белка в 5 компоненте = %.2f\n",
                   recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
        }
    	printf("Аминокислота в продукте = %.2f\n"
    	"Аминокислотный скор = %.2f\n"
    	"Коэффициент рациональности = %.2f\n", akp[i], c[i], aj[i]);
    	printf("Коэффициент сбалансированности = %.2f\n", k[i]);
    	puts("=================================================");
    }
    kras = KRAS(C_SUM(c), min_c);
    bc = BC(kras);
    r = R(aj, akp);
    g = G(akp, min_c, fao_voz2007);
    s = S(k);
    k_general = K_general(s, bc, r);
    printf("Коэффициент разбалансированности = %.2f\n", kras);
    printf("Биологическая ценность = %.2f\n", bc);
    printf("Коэффициент рациональности аминокислотного состава = %.2f\n", r);
    printf("Показатель сопоставимой избыточности = %.2f\n", g);
    printf("Индекс сбалансированности = %.2f\n", s);
    printf("Коэффициент сбалансированности белковой составляющей = %.2f\n\n", k_general);

    puts("=================================================");
    /*ВЫВОД ЛИПИДОВ*/
    for(i = 0; i < LIP_PROP; i++)
    {
        if(comp_num == 1 && i == 0)
            printf("Насыщенных жирных кислот = %.2f\n", recount_lip[0][i]);
        if(comp_num == 1 && i == 1)
            printf("Мононасыщенных жирных кислот = %.2f\n", recount_lip[0][i]);
        if(comp_num == 1 && i == 2)
            printf("Полиненасыщенных жирных кислот = %.2f\n", recount_lip[0][i]);
        if(comp_num == 1 && i == 3)
            printf("Омега-3 = %.2f\n", recount_lip[0][i]);
        if(comp_num == 1 && i == 4)
            printf("Омега-6 = %.2f\n", recount_lip[0][i]);

        if(comp_num == 2)
        {
            if(i == 0)
                printf("Насыщенных жирных кислот в 1 компоненте = %.2f\n"
                       "Насыщенных жирных кислот в 2 компоненте = %.2f\n",
                        recount[0][i], recount[1][i]);
            if(i == 1)
                printf("Мононасыщенных жирных кислот в 1 компоненте = %.2f\n"
                       "Мононасыщенных жирных кислот в 2 компоненте = %.2f\n",
                        recount[0][i], recount[1][i]);
            if(i == 2)
                printf("Полиненасыщенных жирных кислот в 1 компоненте = %.2f\n"
                       "Полиненасыщенных жирных кислот в 2 компоненте = %.2f\n",
                        recount[0][i], recount[1][i]);
            if(i == 3)
                printf("Омега-3 в 1 компоненте = %.2f\n"
                       "Омега-3 в 2 компоненте = %.2f\n",
                        recount[0][i], recount[1][i]);
            if(i == 4)
                printf("Омега-6 в 1 компоненте = %.2f\n"
                       "Омега-6 в 2 компоненте = %.2f\n",
                        recount[0][i], recount[1][i]);
        }
        if(comp_num == 3)
        {
        	if(i == 0)
            printf("Насыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 3 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
        	if(i == 1)
			printf("Мононасыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 3 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
			if(i == 2)
			printf("Полиненасыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 3 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
			if(i == 3)
			printf("Омега-3 в 1 компоненте = %.2f\n"
                   "Омега-3 в 2 компоненте = %.2f\n"
                   "Омега-3 в 3 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
			if(i == 4)
			printf("Омега-6 в 1 компоненте = %.2f\n"
                   "Омега-6 в 2 компоненте = %.2f\n"
                   "Омега-6 в 3 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
        }
        if(comp_num == 4)
        {
            if(i == 0)
            printf("Насыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 3 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 4 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
        	if(i == 1)
			printf("Мононасыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 3 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 4 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
			if(i == 2)
			printf("Полиненасыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 3 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 4 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
			if(i == 3)
			printf("Омега-3 в 1 компоненте = %.2f\n"
                   "Омега-3 в 2 компоненте = %.2f\n"
                   "Омега-3 в 3 компоненте = %.2f\n"
                   "Омега-3 в 4 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
			if(i == 4)
			printf("Омега-6 в 1 компоненте = %.2f\n"
                   "Омега-6 в 2 компоненте = %.2f\n"
                   "Омега-6 в 3 компоненте = %.2f\n"
                   "Омега-6 в 4 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
        }
       if(comp_num == 5)
        {
            if(i == 0)
            printf("Насыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 3 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 4 компоненте = %.2f\n"
                   "Насыщенных жирных кислот в 5 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
        	if(i == 1)
			printf("Мононасыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 3 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 4 компоненте = %.2f\n"
                   "Мононасыщенных жирных кислот в 5 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
			if(i == 2)
			printf("Полиненасыщенных жирных кислот в 1 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 2 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 3 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 4 компоненте = %.2f\n"
                   "Полиненасыщенных жирных кислот в 5 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
			if(i == 3)
			printf("Омега-3 в 1 компоненте = %.2f\n"
                   "Омега-3 в 2 компоненте = %.2f\n"
                   "Омега-3 в 3 компоненте = %.2f\n"
                   "Омега-3 в 4 компоненте = %.2f\n"
                   "Омега-3 в 5 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
			if(i == 4)
			printf("Омега-6 в 1 компоненте = %.2f\n"
                   "Омега-6 в 2 компоненте = %.2f\n"
                   "Омега-6 в 3 компоненте = %.2f\n"
                   "Омега-6 в 4 компоненте = %.2f\n"
                   "Омега-6 в 5 компоненте = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
        }
    	printf("Жирных кислот в 100г продукта = %.2f\n", ratio_calc[i]);
    	puts("=================================================");
    }

    printf("Коэффициент жирнокислотного соответствия (при i=3) = %.2f\n", result1);
    printf("Коэффициент жирнокислотного соответствия (при i=5) = %.2f\n", result2);
    // добавить что то там от 0 до 1
	puts("=================================================");
    printf("Хотите повторить расчёты? Введите 0 если нет и любое другое число если да:\n>");
    scanf(" %d", &ans);
    fflush(stdin);
    }
    while(ans);
    return 0;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/*БЕЛКИ*/
// переменная sign даёт понять функции, какой именно компонент требуется
// пересчитать во время её текущего вызова
double Recount(int comp_num, int sign, double * protein, double comp)
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
double AKP(int comp_num, double * prop, double recount[][AMI], int i)
{
    if(comp_num == 1)
        return recount[0][i];
    if(comp_num == 2)
        return (prop[0] * recount[0][i]) + (prop[1] * recount[1][i]);
    if(comp_num == 3)
        return (prop[0] * recount[0][i]) + (prop[1] * recount[1][i]) + (prop[2] * recount[2][i]);
    if(comp_num == 4)
        return (prop[0] * recount[0][i]) + (prop[1] * recount[1][i]) + (prop[2] * recount[2][i]) + (prop[3] * recount[3][i]);
    if(comp_num == 5)
        return (prop[0] * recount[0][i]) + (prop[1] * recount[1][i]) + (prop[2] * recount[2][i]) + (prop[3] * recount[3][i])
         + (prop[4] * recount[4][i]);
}
//======================================================================================
double Aminoacidskor(double akp, double fao_voz2007)
{
	return akp / fao_voz2007 * 100;
}
//======================================================================================
double Koef_Ration(double min_c, double c)
{
	return min_c / c;
}
//======================================================================================
double C_SUM(double * c)
{
    double sum = 0.0;
    for(int i = 0; i < AMI; i++)
    {
        sum += *c;
        c++;
    }
    return sum;
}
//======================================================================================
double KRAS(double c_sum, double min_c)
{
    return (c_sum - (min_c * AMI)) / AMI;
}
//======================================================================================
double BC(double kras)
{
    return 100 - kras;
}
//======================================================================================
double R(double * aj, double * akp)
{
	double numerator = 0.0; // числитель в формуле
	double denominator = 0.0; // знаменатель в формуле
	for(int i = 0; i < AMI; i++)
	{
		numerator += aj[i] * akp[i];
		denominator += akp[i];
	}
	return numerator / denominator;
}
//======================================================================================
double G(double * akp, double min_c, const double * fao_voz2007)
{
	double numerator = 0.0;
	for(int i = 0; i < AMI; i++)
	{
		numerator += akp[i] - min_c / 100 * fao_voz2007[i];
	}
	return numerator / (min_c / 100);
}
//======================================================================================
double K(double fao_voz2007, double akp)
{
    if(fao_voz2007 <= akp)
        return fao_voz2007 / akp;
    else
        return akp / fao_voz2007;
}
//======================================================================================
double S(double * k)
{
    double temp = 1.0;
    for(int i = 0; i < AMI; i++)
    {
        temp *= k[i];
    }
    return pow(temp, 1.0/9.0);
}
//======================================================================================
double K_general(double s, double bc, double r)
{
    double temp = (s / 1) * (bc / 100) * (r / 1);
    return pow(temp, 1.0/3.0);
}
//======================================================================================
//======================================================================================
//======================================================================================
/*ЛИПИДЫ*/
double Recount_lip(int comp_num, int sign, double * lipids, double ultimate)
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
double Ratio_calc(int comp_num, double * prop, double recount_lip[][LIP_PROP], int i)
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
double Lip_balance_ratio(double fao_voz2008, double recount_lip)
{
	if(recount_lip <= fao_voz2008)
		return recount_lip / fao_voz2008;
	else
		return fao_voz2008 / recount_lip;
}
//======================================================================================
void Fattyacid_compliance(double * Lip_balance_ratio, double* result1, double* result2)
{
	double multi1, multi2 = 1.0;
	for(int i = 0; i < 3; i++)
        multi1 *= Lip_balance_ratio[i];
    *result1 = pow(multi1, 1.0 / 3.0);
	for(int i = 0; i < LIP_PROP; i++)
        multi2 *= Lip_balance_ratio[i];
    *result2 = pow(multi2, 1.0 / 5.0);
}
