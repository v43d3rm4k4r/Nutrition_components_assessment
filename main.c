/*
Daniil Kuprianov, autumn 2019, Saint-Petersburg.

The program evaluates the protein and lipid components
        food products.
            */
#include <stdio.h>
#include <windows.h>
#include <math.h>
#include <stdlib.h>
#include "declarations.h"

#define AMI 9 // количество аминокислот
#define LIP_PROP 5 // липидные хар-ки (аналог аминокислот)
#define MAX_COMP 5 // максимальное количество компонентов

// Данная версия программы будет запускаться только на флешке
// и предназначена для тестовой демонстрации
// #define FLASH_TEST

int main(int argc, char* argv[])
{
    SetConsoleCP(1251);SetConsoleOutputCP(1251);
    SetConsoleTitle("NUTRITION COMPONENTS ASSESSMENT");

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
        double comp[MAX_COMP][AMI] = {{0.0}}; // для аминокислот
        /*ЛИПИДЫ*/
        const double fao_voz2008[MAX_COMP] = {33.33, 33.33, 33.33, 6.67, 26.67};
        double lipids[MAX_COMP] = {0.0}; // кол-во липидов на 100г каждого компонента
        /*Здесь один массив отвечает за ненасыщенные жирные кислоты,
        мононенасыщенные жирные кислоты, полиненасыщенные жирные кислоты, омега3 и омега6.
        Это сделано для того чтобы передать один двумерный массив в функцию вместо 5 обычных.
        [0] - насыщенные жирные кислоты
        [1] - мононенасыщенные жирные кислоты
        [2] - полиненасыщенные жирные кислоты
        [3] - омега6
        [4] - омега3*/
        double ultimate[MAX_COMP][LIP_PROP] = {{0.0}};

        /*ДЛЯ РЕЗУЛЬТАТОВ ВЫЧИСЛЕНИЙ*/
        /*БЕЛКИ*/
        double akp[AMI] = {0.0};
        double aminoacidskor[AMI] = {0.0};
        double koef_ration[AMI] = {0.0};
        double fatty_acid_per_100g[AMI] = {0.0};
        double prop[MAX_COMP] = {0.0}; // для пропорций, если компонентов > 1
        double recount[MAX_COMP][AMI] = {{0.0}}; // для пересчёта аминокилот
        /*коэффициент разбалансированности, биологическая ценность,
        коэффициент рациональности аминокислотного состава,
        индекс сбалансированности, показатель сопоставимой избыточности и
        коэффициент сбалансированности всей белковой составляющей(вычисляются 1 раз в конце)*/
        double kras = 0.0;
        double biological_value = 0.0; // bc
        double amino_acid_comp_ratio_coef = 0.0; // r
        double comparable_redundancy_ratio = 0.0; // g
        double balance_index = 0.0; // s
        double k_general = 0.0;
        int i; // счётчик, номер строки с аминокислотой

        /*ЛИПИДЫ*/
        double recount_lip[MAX_COMP][LIP_PROP] = {{0.0}};
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
        [3] - омега6
        [4] - омега3*/
        for(i = 0; i < LIP_PROP; i++)
        {
            if(comp_num > 1 && i == 0)
                puts("Введите кол-во насыщенных жирных кислот на 100г каждого компонента:");
            if(comp_num > 1 && i == 1)
                puts("Введите кол-во мононенасыщенных жирных кислот на 100г каждого компонента:");
            if(comp_num > 1 && i == 2)
                puts("Введите кол-во полиненасыщенных жирных кислот на 100г каждого компонента:");
            if(comp_num > 1 && i == 3)
                puts("Введите омега-6 на 100г каждого компонента:");
            if(comp_num > 1 && i == 4)
                puts("Введите омега-3 на 100г каждого компонента:");
            if(comp_num == 1)
            {
                if(i == 0)
                    puts("Введите кол-во насыщенных жирных кислот на 100г компонента:");
                if(i == 1)
                    puts("Введите кол-во мононенасыщенных жирных кислот на 100г компонента:");
                if(i == 2)
                    puts("Введите кол-во полиненасыщенных жирных кислот на 100г компонента:");
                if(i == 3)
                    puts("Введите кол-во омега-6 на 100г компонента:");
                if(i == 4)
                    puts("Введите кол-во омега-3 на 100г компонента:");

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
            aminoacidskor[i] = Aminoacidskor(akp[i], fao_voz2007[i]);
            fatty_acid_per_100g[i] = Fatty_Acid_Per_100g(fao_voz2007[i], akp[i]);
        }

        // в конце находим минимальный 'с'(aminoacidskor) из всех, считаем aj для каждой аминокислоты,
        // коэффициент рациональности аминокислотного состава и показатель сопоставимой избыточности
        double min_c = 999.0;
        for(i = 0; i < AMI; i++)
        {
            if(aminoacidskor[i] < min_c)
                min_c = aminoacidskor[i];
        }
        for(i = 0; i < AMI; i++)
        {
            koef_ration[i] = Koef_Ration(min_c, aminoacidskor[i]);
        }
        //======================================================================================
        /*ВЫЧИСЛЕНИЯ ЛИПИДОВ*/
        // здесь i это номер "свойства" липидов
        for(i = 0; i < LIP_PROP; i++)
        {
            if(comp_num == 1)
                recount_lip[0][i] = Recount_Lip(comp_num, 1, lipids, ultimate[0][i]);
            if(comp_num == 2)
            {
                recount_lip[0][i] = Recount_Lip(comp_num, 1, lipids, ultimate[0][i]);
                recount_lip[1][i] = Recount_Lip(comp_num, 2, lipids, ultimate[1][i]);
            }
            if(comp_num == 3)
            {
                recount_lip[0][i] = Recount_Lip(comp_num, 1, lipids, ultimate[0][i]);
                recount_lip[1][i] = Recount_Lip(comp_num, 2, lipids, ultimate[1][i]);
                recount_lip[2][i] = Recount_Lip(comp_num, 3, lipids, ultimate[2][i]);
            }
            if(comp_num == 4)
            {
                recount_lip[0][i] = Recount_Lip(comp_num, 1, lipids, ultimate[0][i]);
                recount_lip[1][i] = Recount_Lip(comp_num, 2, lipids, ultimate[1][i]);
                recount_lip[2][i] = Recount_Lip(comp_num, 3, lipids, ultimate[2][i]);
                recount_lip[3][i] = Recount_Lip(comp_num, 4, lipids, ultimate[3][i]);
            }
            if(comp_num == 5)
            {
                recount_lip[0][i] = Recount_Lip(comp_num, 1, lipids, ultimate[0][i]);
                recount_lip[1][i] = Recount_Lip(comp_num, 2, lipids, ultimate[1][i]);
                recount_lip[2][i] = Recount_Lip(comp_num, 3, lipids, ultimate[2][i]);
                recount_lip[3][i] = Recount_Lip(comp_num, 4, lipids, ultimate[3][i]);
                recount_lip[4][i] = Recount_Lip(comp_num, 5, lipids, ultimate[4][i]);
            }
            ratio_calc[i] = Ratio_Calc(comp_num, prop, recount_lip, i);
            lip_balance_ratio[i] = Lip_Balance_Ratio(fao_voz2008[i], ratio_calc[i]);
            Fattyacid_Compliance(lip_balance_ratio, &result1, &result2);
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
                   "Коэффициент рациональности = %.2f\n", akp[i], aminoacidskor[i], koef_ration[i]);
            printf("Коэффициент сбалансированности = %.2f\n", fatty_acid_per_100g[i]);
            puts("=================================================");
        }
        kras = KRAS(Aminoacidskor_Sum(aminoacidskor), min_c);
        biological_value = Biological_Value(kras);
        amino_acid_comp_ratio_coef = Amino_Acid_Comp_Ratio_Coef(koef_ration, akp);
        comparable_redundancy_ratio = Comparable_Redundancy_Ratio(akp, min_c, fao_voz2007);
        balance_index = Balance_Index(fatty_acid_per_100g);
        k_general = Balance_Index_General(balance_index, biological_value, amino_acid_comp_ratio_coef);
        puts("==========ОЦЕНКА БЕЛКОВОЙ СОСТАВЛЯЮЩЕЙ:==========\n");
        printf("Коэффициент разбалансированности = %.2f\n", kras);
        printf("Биологическая ценность = %.2f\n", biological_value);
        printf("Коэффициент рациональности аминокислотного состава = %.2f\n", amino_acid_comp_ratio_coef);
        printf("Показатель сопоставимой избыточности = %.2f\n", comparable_redundancy_ratio);
        printf("Индекс сбалансированности = %.2f\n", balance_index);
        printf("Коэффициент сбалансированности белковой составляющей = %.2f\n", k_general);

        if(k_general >= 0.0 && k_general <= 0.2)
            printf("Шкала желательности Харрингтона - очень низкая\n");
        if(k_general > 0.2 && k_general <= 0.37)
            printf("Шкала желательности Харрингтона - низкая\n");
        if(k_general > 0.37 && k_general <= 0.64)
            printf("Шкала желательности Харрингтона - средняя\n");
        if(k_general > 0.64 && k_general <= 0.8)
            printf("Шкала желательности Харрингтона - высокая\n");
        if(k_general > 0.8 && k_general <= 1.0)
            printf("Шкала желательности Харрингтона - очень высокая\n");

        puts("=================================================");
        /*ВЫВОД ЛИПИДОВ*/
        for(i = 0; i < LIP_PROP; i++)
        {
            if(comp_num == 1 && i == 0)
                printf("Насыщенных жирных кислот = %.2f\n", recount_lip[0][i]);
            if(comp_num == 1 && i == 1)
                printf("Мононенасыщенных жирных кислот = %.2f\n", recount_lip[0][i]);
            if(comp_num == 1 && i == 2)
                printf("Полиненасыщенных жирных кислот = %.2f\n", recount_lip[0][i]);
            if(comp_num == 1 && i == 3)
                printf("Омега-6 = %.2f\n", recount_lip[0][i]);
            if(comp_num == 1 && i == 4)
                printf("Омега-3 = %.2f\n", recount_lip[0][i]);

            if(comp_num == 2)
            {
                if(i == 0)
                    printf("Насыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 2 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i]);
                if(i == 1)
                    printf("Мононенасыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 2 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i]);
                if(i == 2)
                    printf("Полиненасыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 2 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i]);
                if(i == 3)
                    printf("Омега-6 в 1 компоненте = %.2f\n"
                           "Омега-6 в 2 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i]);
                if(i == 4)
                    printf("Омега-3 в 1 компоненте = %.2f\n"
                           "Омега-3 в 2 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i]);
            }
            if(comp_num == 3)
            {
                if(i == 0)
                    printf("Насыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 3 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i]);
                if(i == 1)
                    printf("Мононенасыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 3 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i]);
                if(i == 2)
                    printf("Полиненасыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 3 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i]);
                if(i == 3)
                    printf("Омега-6 в 1 компоненте = %.2f\n"
                           "Омега-6 в 2 компоненте = %.2f\n"
                           "Омега-6 в 3 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i]);
                if(i == 4)
                    printf("Омега-3 в 1 компоненте = %.2f\n"
                           "Омега-3 в 2 компоненте = %.2f\n"
                           "Омега-3 в 3 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i]);
            }
            if(comp_num == 4)
            {
                if(i == 0)
                    printf("Насыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 3 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 4 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i]);
                if(i == 1)
                    printf("Мононенасыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 3 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 4 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i]);
                if(i == 2)
                    printf("Полиненасыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 3 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 4 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i]);
                if(i == 3)
                    printf("Омега-6 в 1 компоненте = %.2f\n"
                           "Омега-6 в 2 компоненте = %.2f\n"
                           "Омега-6 в 3 компоненте = %.2f\n"
                           "Омега-6 в 4 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i]);
                if(i == 4)
                    printf("Омега-3 в 1 компоненте = %.2f\n"
                           "Омега-3 в 2 компоненте = %.2f\n"
                           "Омега-3 в 3 компоненте = %.2f\n"
                           "Омега-3 в 4 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i]);
            }
            if(comp_num == 5)
            {
                if(i == 0)
                    printf("Насыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 3 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 4 компоненте = %.2f\n"
                           "Насыщенных жирных кислот в 5 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i], recount_lip[4][i]);
                if(i == 1)
                    printf("Мононенасыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 3 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 4 компоненте = %.2f\n"
                           "Мононенасыщенных жирных кислот в 5 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i], recount_lip[4][i]);
                if(i == 2)
                    printf("Полиненасыщенных жирных кислот в 1 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 2 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 3 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 4 компоненте = %.2f\n"
                           "Полиненасыщенных жирных кислот в 5 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i], recount_lip[4][i]);
                if(i == 3)
                    printf("Омега-6 в 1 компоненте = %.2f\n"
                           "Омега-6 в 2 компоненте = %.2f\n"
                           "Омега-6 в 3 компоненте = %.2f\n"
                           "Омега-6 в 4 компоненте = %.2f\n"
                           "Омега-6 в 5 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i], recount_lip[4][i]);
                if(i == 4)
                    printf("Омега-3 в 1 компоненте = %.2f\n"
                           "Омега-3 в 2 компоненте = %.2f\n"
                           "Омега-3 в 3 компоненте = %.2f\n"
                           "Омега-3 в 4 компоненте = %.2f\n"
                           "Омега-3 в 5 компоненте = %.2f\n",
                           recount_lip[0][i], recount_lip[1][i], recount_lip[2][i], recount_lip[3][i], recount_lip[4][i]);
            }
            printf("Жирных кислот в 100г продукта = %.2f\n", ratio_calc[i]);
            puts("=================================================");
        }
        puts("==========ОЦЕНКА ЛИПИДНОЙ СОСТАВЛЯЮЩЕЙ:==========\n");
        printf("Коэффициент жирнокислотного соответствия (при i=3) = %.2f\n", result1);
        printf("Коэффициент жирнокислотного соответствия (при i=5) = %.2f\n", result2);

        if(result1 >= 0.0 && result1 <= 0.2)
            printf("Шкала желательности Харрингтона (при i=3) - очень низкая\n");
        if(result1 > 0.2 &&result1 <= 0.37)
            printf("Шкала желательности Харрингтона (при i=3) - низкая\n");
        if(result1 > 0.37 && result1 <= 0.64)
            printf("Шкала желательности Харрингтона (при i=3) - средняя\n");
        if(result1 > 0.64 && result1 <= 0.8)
            printf("Шкала желательности Харрингтона (при i=3) - высокая\n");
        if(result1 > 0.8 && result1 <= 1.0)
            printf("Шкала желательности Харрингтона (при i=3) - очень высокая\n");

        if(result2 >= 0.0 && result2 <= 0.2)
            printf("Шкала желательности Харрингтона (при i=5) - очень низкая\n");
        if(result2 > 0.2 &&result2 <= 0.37)
            printf("Шкала желательности Харрингтона (при i=5) - низкая\n");
        if(result2 > 0.37 && result2 <= 0.64)
            printf("Шкала желательности Харрингтона (при i=5) - средняя\n");
        if(result2 > 0.64 && result2 <= 0.8)
            printf("Шкала желательности Харрингтона (при i=5) - высокая\n");
        if(result2 > 0.8 && result2 <= 1.0)
            printf("Шкала желательности Харрингтона (при i=5) - очень высокая\n");

        puts("=================================================");
        printf("Хотите повторить расчёты? Введите 0 если нет и любое другое число если да:\n>");
        scanf(" %d", &ans);
        fflush(stdin);
    }
    while(ans);
    return 0;
}
