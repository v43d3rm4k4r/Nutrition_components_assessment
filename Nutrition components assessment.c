/*������������ ��������� �.�., ����� 2019 ����.
 ��������� ���������� ������ ��������  � �������� ������������
 ��������� �������.*/
 /*
 1.������ ���-�� �����������
 2.������ ���-�� ����� �� 100� ������� ���������� (����� ����� ���-�� �����������) (n � �������)(�� 1 �� 5, �������� �� 0 �� 100)
 3.������ ��������� �����������
 4.������ ������������ 9 ��� ��� ������� ���������� (X���� � �������)
 5.������������� ������������ � �������� �� �����
 4.��������� ���, �, �����, ����, ��, R, G
*/
#include <stdio.h>
#include <windows.h>
#include <math.h>
#include <stdlib.h>

#define AMI 9 // ���������� �����������
#define LIP_PROP 5 // �������� ���-�� (������ �����������)
#define MAX_COMP 5 // ������������ ���������� �����������

// ������ ������ ��������� ����� ����������� ������ �� ������
// � ������������� ��� �������� ������������
// #define FLASH_TEST
//======================================================================================
/*�����*/
// �������� ����������� �� 100� ���������� �� 100� �����
double Recount(int comp_num, int sign, double * protein, double comp);
// ������������ � ��������
double AKP(int comp_num, double * prop, double recount[][AMI], int i);
// �������������� ����
double Aminoacidskor(double akp, const double fao_voz2007);
// ����������� ��������������
double Koef_Ration(double min_c, double c);
// ����� �������������� ������ ��� ���� �����������
double C_SUM(double * c);
// ����������� ��������������������
double KRAS(double c_sum, double min_c);
// ������������� ��������
double BC(double kras);
// ����������� �������������� ��������������� �������
double R(double * aj, double * akp);
// ���������� ������������ ������������
double G(double * akp, double min_c, const double * fao_voz2007);
// ������� � ������ ��������� (������ ������� �� 100� ��������, ������ ���)
double K(const double fao_voz2007, double akp);
// ����������� ������������������
double S(double * k);
// ����������� ������������������ ���� �������� ������������
double K_general(double s, double bc, double r);
//======================================================================================
/*������*/
// �������� �� 100� ���� ���� �������� ������� ������
double Recount_lip(int comp_num, int sign, double * lipids, double ultimate);
// ����������� �������������� ������ �������
double Ratio_calc(int comp_num, double * prop, double recount_lip[][LIP_PROP], int i);
// ����������� ������������������
double Lip_balance_ratio(double fao_voz2008, double recount_lip);
// ����������� ��������������� ������������
void Fattyacid_compliance(double * Lip_balance_ratio, double* result1, double* result2);

int main(void)
{
	SetConsoleCP(1251);SetConsoleOutputCP(1251);
	#ifdef FLASH_TEST
	int is_flash = GetDriveType(NULL); // NULL ������������� ��� ����������� ������� ��������
	if(is_flash != 2)
    {
        printf("������ �������!");
        exit(EXIT_FAILURE);
    }
    #endif // FLASH_TEST
    int ans = 0;
    system("color 70");
    puts("�� ���� �������� - vaedermakar@mail.ru.\n");
    puts("*��������� ���������� ��� ��� 2007 ��� ���������� �������� ������������\n"
        "��� ��� 2008 ��� ���������� �������� ������������.");
    do{
    /*��� �������� ������*/
    /*�����*/
    int comp_num = 0; // ���-�� �����������
    /*�����*/
    const double fao_voz2007[AMI] = {3.9, 1.5, 3.0, 5.9, 4.5, 2.2, 2.3, 0.6, 3.8};
    double protein[MAX_COMP] = {0.0}; // ���-�� ����� �� 100� ������� ����������
	double comp[MAX_COMP][AMI] = {0.0}; // ��� �����������
	/*������*/
    const double fao_voz2008[MAX_COMP] = {33.33, 33.33, 33.33, 6.67, 26.67};
    double lipids[MAX_COMP] = {0.0}; // ���-�� ������� �� 100� ������� ����������
    /*����� ���� ������ �������� �� ������������ ������ �������,
    ���������������� ������ �������, ���������������� ������ �������, �����3 � �����6.
    ��� ������� ��� ���� ����� �������� ���� ��������� ������ � ������� ������ 5 �������.
    [0] - ���������� ������ �������
    [1] - ���������������� ������ �������
    [2] - ���������������� ������ �������
    [3] - �����3
    [4] - �����6*/
    double ultimate[MAX_COMP][LIP_PROP] = {0.0};

	/*��� ����������� ����������*/
	/*�����*/
	double akp[AMI] = {0.0};
	double c[AMI] = {0.0};
	double aj[AMI] = {0.0};
	double k[AMI] = {0.0};
	double prop[MAX_COMP] = {0.0}; // ��� ���������, ���� ����������� > 1
	double recount[MAX_COMP][AMI] = {0.0}; // ��� ��������� ����������
    /*����������� ��������������������, ������������� ��������,
    ����������� �������������� ��������������� �������,
    ������ ������������������, ���������� ������������ ������������ �
    ����������� ������������������ ���� �������� ������������(����������� 1 ��� � �����)*/
    double kras = 0.0;
    double bc = 0.0;
    double r = 0.0;
    double g = 0.0;
    double s = 0.0;
    double k_general = 0.0;
    int i; // �������, ����� ������ � �������������

    /*������*/
    double recount_lip[MAX_COMP][LIP_PROP] = {0.0};
    double ratio_calc[LIP_PROP] = {0.0};
	double lip_balance_ratio[LIP_PROP] = {0.0};
	double result1 = 0.0;
	double result2 = 0.0;
    //======================================================================================
    do{
    puts("������� ���������� ����������� �������� (�� 5):");
    while(scanf(" %d", &comp_num) != 1)
    {
        {
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    putchar('\n');
    // ���� ����������� 0 ��� ������ ���� ������ 5
    if(comp_num <= 0 || comp_num > MAX_COMP)
        printf("����������� �� ����� ���� %d, ���������� ��� ���.\n", comp_num);
    }while(comp_num <= 0 || comp_num > MAX_COMP);
    //======================================================================================
    do{
    if(comp_num == 1)
    {
        prop[0] = 1.0;
    }
    if(comp_num > 1)
    {
        puts("������� ��������� �����������(����� ������� - 0.20 0.60 0.20):");
        if(comp_num == 2)
        {
            while(scanf(" %lf %lf", &prop[0], &prop[1]) != 2)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
        }
        if(comp_num == 3)
        {
            while(scanf(" %lf %lf %lf", &prop[0], &prop[1], &prop[2]) != 3)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
        }
        if(comp_num == 4)
        {
            while(scanf(" %lf %lf %lf %lf", &prop[0], &prop[1], &prop[2], &prop[3]) != 4)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
        }
        if(comp_num == 5)
        {
            while(scanf(" %lf %lf %lf %lf %lf", &prop[0], &prop[1], &prop[2], &prop[3], &prop[4]) != 5)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
        }
        putchar('\n');
    }
    // ������ �������� ���������� � �������
    // ���� ����������� ������ 1 � �� ������ 5, ����� ��������� ������ ��������� 1.0(100%)
    if(fabs((prop[0]+prop[1]+prop[2]+prop[3]+prop[4])-(1.0)) >= 0.000000000001)
        puts("������ �����. ����� ��������� ����������� ������ ��������� 1.0.");
	}while(fabs((prop[0]+prop[1]+prop[2]+prop[3]+prop[4])-(1.0)) >= 0.000000000001);
    //======================================================================================
    if(comp_num > 1)
        puts("������� ���-�� ����� �� 100� ������� ����������:");
    if(comp_num == 1)
    {
        puts("������� ���������� ����� �� 100� ����������:");
        while(scanf(" %lf", &protein[0]) != 1)
        {
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    if(comp_num == 2)
    {
        while(scanf(" %lf %lf", &protein[0], &protein[1]) != 2
            || protein[0] < 0.0 || protein[0] > 100.0 || protein[1] < 0.0 || protein[1] > 100.0)
        {
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    if(comp_num == 3)
    {
        while(scanf(" %lf %lf %lf", &protein[0], &protein[1], &protein[2]) != 3
            || protein[0] < 0.0 || protein[0] > 100.0 || protein[1] < 0.0 || protein[1] > 100.0
            || protein[2] < 0.0 || protein[2] > 100.0)
        {
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    if(comp_num == 4)
    {
        while(scanf(" %lf %lf %lf %lf", &protein[0], &protein[1], &protein[2], &protein[3]) != 4
            || protein[0] < 0.0 || protein[0] > 100.0 || protein[1] < 0.0 || protein[1] > 100.0
            || protein[2] < 0.0 || protein[2] > 100.0 || protein[3] < 0.0 || protein[3] > 100.0)
        {
            puts("������������ ����.");
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
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    putchar('\n');
    //======================================================================================
    if(comp_num > 1)
        puts("������� ������������ � ������� ��� ������� ����������.");
    if(comp_num == 1)
        puts("������� ������������ � �������.");
    for(i = 0; i < AMI; i++)
    {
        if(i == 0 && comp_num > 1)
            puts("������� ����� ��� ������� ����������:");
        if(i == 1 && comp_num > 1)
            puts("������� �������� ��� ������� ����������:");
        if(i == 2 && comp_num > 1)
           puts("������� ��������� ��� ������� ����������:");
        if(i == 3 && comp_num > 1)
           puts("������� ������ ��� ������� ����������:");
        if(i == 4 && comp_num > 1)
            puts("������� ����� ��� ������� ����������:");
        if(i == 5 && comp_num > 1)
            puts("������� ��������+������� ��� ������� ����������:");
        if(i == 6 && comp_num > 1)
            puts("������� ������� ��� ������� ����������:");
        if(i == 7 && comp_num > 1)
            puts("������� ��������� ��� ������� ����������:");
        if(i == 8 && comp_num > 1)
           puts("������� �����������+������� ��� ������� ����������:");

        if(i == 0 && comp_num == 1)
            puts("������� �����:");
        if(i == 1 && comp_num == 1)
            puts("������� ��������:");
        if(i == 2 && comp_num == 1)
            puts("������� ���������:");
        if(i == 3 && comp_num == 1)
            puts("������� ������:");
        if(i == 4 && comp_num == 1)
            puts("������� �����:");
        if(i == 5 && comp_num == 1)
            puts("������� ��������+�������:");
        if(i == 6 && comp_num == 1)
            puts("������� �������:");
        if(i == 7 && comp_num == 1)
            puts("������� ���������:");
        if(i == 8 && comp_num == 1)
            puts("������� �����������+�������:");

        if(comp_num == 1)
        {
            while(scanf(" %lf", &comp[0][i]) != 1 || comp[0][i] <= 0.0)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 2)
        {
            while(scanf(" %lf %lf", &comp[0][i], &comp[1][i]) != 2
                  || comp[0][i] <= 0.0 || comp[1][i] <= 0.0)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 3)
        {
            while(scanf(" %lf %lf %lf", &comp[0][i], &comp[1][i], &comp[2][i]) != 3
                  || comp[0][i] <= 0.0 || comp[1][i] <= 0.0 || comp[2][i] <= 0.0)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
            putchar('\n');
        }
    	if(comp_num == 4)
        {
            while(scanf(" %lf %lf %lf %lf", &comp[0][i], &comp[1][i], &comp[2][i], &comp[3][i]) != 4
                  || comp[0][i] <= 0.0 || comp[1][i] <= 0.0 || comp[2][i] <= 0.0 || comp[3][i] <= 0.0)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 5)
        {
            while(scanf(" %lf %lf %lf %lf %lf", &comp[0][i], &comp[1][i], &comp[2][i], &comp[3][i], &comp[4][i]) != 5
                  || comp[0][i] <= 0.0 || comp[1][i] <= 0.0 || comp[2][i] <= 0.0 || comp[3][i] <= 0.0 || comp[4][i] <= 0.0)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
            putchar('\n');
        }
    }
    //======================================================================================
    // ���� ����� �� 100� ������� ���������� ��� ���������
    if(comp_num > 1)
        puts("������� ���-�� ������� �� 100� ������� ����������:");
    if(comp_num == 1)
    {
        puts("������� ���-�� ������� �� 100� ����������:");
        while(scanf(" %lf", &lipids[0]) != 1 || lipids[0] < 0.0 || lipids[0] > 100.0)
        {
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    if(comp_num == 2)
    {
        while(scanf(" %lf %lf", &lipids[0], &lipids[1]) != 2
            || lipids[0] < 0.0 || lipids[0] > 100.0 || lipids[1] < 0.0 || lipids[1] > 100.0)
        {
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    if(comp_num == 3)
    {
        while(scanf(" %lf %lf %lf", &lipids[0], &lipids[1], &lipids[2]) != 3
            || lipids[0] < 0.0 || lipids[0] > 100.0 || lipids[1] < 0.0 || lipids[1] > 100.0
            || lipids[2] < 0.0 || lipids[2] > 100.0)
        {
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    if(comp_num == 4)
    {
        while(scanf(" %lf %lf %lf %lf", &lipids[0], &lipids[1], &lipids[2], &lipids[3]) != 4
            || lipids[0] < 0.0 || lipids[0] > 100.0 || lipids[1] < 0.0 || lipids[1] > 100.0
            || lipids[2] < 0.0 || lipids[2] > 100.0 || lipids[3] < 0.0 || lipids[3] > 100.0)
        {
            puts("������������ ����.");
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
            puts("������������ ����.");
            fflush(stdin);
        }
    }
    putchar('\n');
    //======================================================================================
   /*����:
    [0] - ���������� ������ �������
    [1] - ���������������� ������ �������
    [2] - ���������������� ������ �������
    [3] - �����3
    [4] - �����6*/
    for(i = 0; i < LIP_PROP; i++)
    {
        if(comp_num > 1 && i == 0)
            puts("������� ���-�� ���������� ������ ������ �� 100� ������� ����������:");
        if(comp_num > 1 && i == 1)
            puts("������� ���-�� ���������������� ������ ������ �� 100� ������� ����������:");
        if(comp_num > 1 && i == 2)
            puts("������� ���-�� ���������������� ������ ������ �� 100� ������� ����������:");
        if(comp_num > 1 && i == 3)
            puts("������� �����-3 �� 100� ������� ����������:");
        if(comp_num > 1 && i == 4)
            puts("������� �����-6 �� 100� ������� ����������:");
        if(comp_num == 1)
        {
            if(i == 0)
                puts("������� ���-�� ���������� ������ ������ �� 100� ����������:");
            if(i == 1)
                puts("������� ���-�� ���������������� ������ ������ �� 100� ����������:");
            if(i == 2)
                puts("������� ���-�� ���������������� ������ ������ �� 100� ����������:");
            if(i == 3)
                puts("������� ���-�� �����-3 �� 100� ����������:");
            if(i == 4)
                puts("������� ���-�� �����-6 �� 100� ����������:");

            while(scanf(" %lf", &ultimate[0][i]) != 1 || ultimate[0][i] < 0.0 || ultimate[0][i] > 100.0)
            {
                puts("������������ ����.");
                fflush(stdin);
            }
            putchar('\n');
        }
        if(comp_num == 2)
        {
            while(scanf(" %lf %lf", &ultimate[0][i], &ultimate[1][i]) != 2
                || ultimate[0][i] < 0.0 || ultimate[0][i] > 100.0 || ultimate[1][i] < 0.0 || ultimate[1][i] > 100.0)
            {
                puts("������������ ����.");
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
                puts("������������ ����.");
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
                puts("������������ ����.");
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
                puts("������������ ����.");
                fflush(stdin);
            }
            putchar('\n');
        }
    }
    //======================================================================================
    /*���������� ������*/
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

    // � ����� ������� ����������� '�' �� ����, ������� aj ��� ������ ������������,
    // ����������� �������������� ��������������� ������� � ���������� ������������ ������������
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
    /*���������� �������*/
    // ����� i ��� ����� "��������" �������
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
    /*����� �����������*/
    puts("=================================================");
    for(i = 0; i < AMI; i++)
    {
        if(i == 0)
            puts("�����:\n");
        if(i == 1)
            puts("��������:\n");
        if(i == 2)
            puts("���������:\n");
        if(i == 3)
            puts("������:\n");
        if(i == 4)
            puts("�����:\n");
        if(i == 5)
            puts("��������+�������:\n");
        if(i == 6)
            puts("�������:\n");
        if(i == 7)
            puts("���������:\n");
        if(i == 8)
            puts("�����������+�������:\n");

    	if(comp_num == 1)
            printf("������������ �� 100� ����� = %.2f\n", recount[0][i]);
        if(comp_num == 2)
        {
            printf("������������ �� 100� ����� � 1 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 2 ���������� = %.2f\n",
                    recount[0][i], recount[1][i]);
        }
        if(comp_num == 3)
        {
            printf("������������ �� 100� ����� � 1 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 2 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 3 ���������� = %.2f\n",
                   recount[0][i], recount[1][i], recount[2][i]);
        }
        if(comp_num == 4)
        {
            printf("������������ �� 100� ����� � 1 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 2 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 3 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 4 ���������� = %.2f\n",
                   recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
        }
        if(comp_num == 5)
        {
            printf("������������ �� 100� ����� � 1 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 2 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 3 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 4 ���������� = %.2f\n"
                   "������������ �� 100� ����� � 5 ���������� = %.2f\n",
                   recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
        }
    	printf("������������ � �������� = %.2f\n"
    	"�������������� ���� = %.2f\n"
    	"����������� �������������� = %.2f\n", akp[i], c[i], aj[i]);
    	printf("����������� ������������������ = %.2f\n", k[i]);
    	puts("=================================================");
    }
    kras = KRAS(C_SUM(c), min_c);
    bc = BC(kras);
    r = R(aj, akp);
    g = G(akp, min_c, fao_voz2007);
    s = S(k);
    k_general = K_general(s, bc, r);
    printf("����������� �������������������� = %.2f\n", kras);
    printf("������������� �������� = %.2f\n", bc);
    printf("����������� �������������� ��������������� ������� = %.2f\n", r);
    printf("���������� ������������ ������������ = %.2f\n", g);
    printf("������ ������������������ = %.2f\n", s);
    printf("����������� ������������������ �������� ������������ = %.2f\n\n", k_general);

    puts("=================================================");
    /*����� �������*/
    for(i = 0; i < LIP_PROP; i++)
    {
        if(comp_num == 1 && i == 0)
            printf("���������� ������ ������ = %.2f\n", recount_lip[0][i]);
        if(comp_num == 1 && i == 1)
            printf("�������������� ������ ������ = %.2f\n", recount_lip[0][i]);
        if(comp_num == 1 && i == 2)
            printf("���������������� ������ ������ = %.2f\n", recount_lip[0][i]);
        if(comp_num == 1 && i == 3)
            printf("�����-3 = %.2f\n", recount_lip[0][i]);
        if(comp_num == 1 && i == 4)
            printf("�����-6 = %.2f\n", recount_lip[0][i]);

        if(comp_num == 2)
        {
            if(i == 0)
                printf("���������� ������ ������ � 1 ���������� = %.2f\n"
                       "���������� ������ ������ � 2 ���������� = %.2f\n",
                        recount[0][i], recount[1][i]);
            if(i == 1)
                printf("�������������� ������ ������ � 1 ���������� = %.2f\n"
                       "�������������� ������ ������ � 2 ���������� = %.2f\n",
                        recount[0][i], recount[1][i]);
            if(i == 2)
                printf("���������������� ������ ������ � 1 ���������� = %.2f\n"
                       "���������������� ������ ������ � 2 ���������� = %.2f\n",
                        recount[0][i], recount[1][i]);
            if(i == 3)
                printf("�����-3 � 1 ���������� = %.2f\n"
                       "�����-3 � 2 ���������� = %.2f\n",
                        recount[0][i], recount[1][i]);
            if(i == 4)
                printf("�����-6 � 1 ���������� = %.2f\n"
                       "�����-6 � 2 ���������� = %.2f\n",
                        recount[0][i], recount[1][i]);
        }
        if(comp_num == 3)
        {
        	if(i == 0)
            printf("���������� ������ ������ � 1 ���������� = %.2f\n"
                   "���������� ������ ������ � 2 ���������� = %.2f\n"
                   "���������� ������ ������ � 3 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
        	if(i == 1)
			printf("�������������� ������ ������ � 1 ���������� = %.2f\n"
                   "�������������� ������ ������ � 2 ���������� = %.2f\n"
                   "�������������� ������ ������ � 3 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
			if(i == 2)
			printf("���������������� ������ ������ � 1 ���������� = %.2f\n"
                   "���������������� ������ ������ � 2 ���������� = %.2f\n"
                   "���������������� ������ ������ � 3 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
			if(i == 3)
			printf("�����-3 � 1 ���������� = %.2f\n"
                   "�����-3 � 2 ���������� = %.2f\n"
                   "�����-3 � 3 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
			if(i == 4)
			printf("�����-6 � 1 ���������� = %.2f\n"
                   "�����-6 � 2 ���������� = %.2f\n"
                   "�����-6 � 3 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i]);
        }
        if(comp_num == 4)
        {
            if(i == 0)
            printf("���������� ������ ������ � 1 ���������� = %.2f\n"
                   "���������� ������ ������ � 2 ���������� = %.2f\n"
                   "���������� ������ ������ � 3 ���������� = %.2f\n"
                   "���������� ������ ������ � 4 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
        	if(i == 1)
			printf("�������������� ������ ������ � 1 ���������� = %.2f\n"
                   "�������������� ������ ������ � 2 ���������� = %.2f\n"
                   "�������������� ������ ������ � 3 ���������� = %.2f\n"
                   "�������������� ������ ������ � 4 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
			if(i == 2)
			printf("���������������� ������ ������ � 1 ���������� = %.2f\n"
                   "���������������� ������ ������ � 2 ���������� = %.2f\n"
                   "���������������� ������ ������ � 3 ���������� = %.2f\n"
                   "���������������� ������ ������ � 4 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
			if(i == 3)
			printf("�����-3 � 1 ���������� = %.2f\n"
                   "�����-3 � 2 ���������� = %.2f\n"
                   "�����-3 � 3 ���������� = %.2f\n"
                   "�����-3 � 4 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
			if(i == 4)
			printf("�����-6 � 1 ���������� = %.2f\n"
                   "�����-6 � 2 ���������� = %.2f\n"
                   "�����-6 � 3 ���������� = %.2f\n"
                   "�����-6 � 4 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i]);
        }
       if(comp_num == 5)
        {
            if(i == 0)
            printf("���������� ������ ������ � 1 ���������� = %.2f\n"
                   "���������� ������ ������ � 2 ���������� = %.2f\n"
                   "���������� ������ ������ � 3 ���������� = %.2f\n"
                   "���������� ������ ������ � 4 ���������� = %.2f\n"
                   "���������� ������ ������ � 5 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
        	if(i == 1)
			printf("�������������� ������ ������ � 1 ���������� = %.2f\n"
                   "�������������� ������ ������ � 2 ���������� = %.2f\n"
                   "�������������� ������ ������ � 3 ���������� = %.2f\n"
                   "�������������� ������ ������ � 4 ���������� = %.2f\n"
                   "�������������� ������ ������ � 5 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
			if(i == 2)
			printf("���������������� ������ ������ � 1 ���������� = %.2f\n"
                   "���������������� ������ ������ � 2 ���������� = %.2f\n"
                   "���������������� ������ ������ � 3 ���������� = %.2f\n"
                   "���������������� ������ ������ � 4 ���������� = %.2f\n"
                   "���������������� ������ ������ � 5 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
			if(i == 3)
			printf("�����-3 � 1 ���������� = %.2f\n"
                   "�����-3 � 2 ���������� = %.2f\n"
                   "�����-3 � 3 ���������� = %.2f\n"
                   "�����-3 � 4 ���������� = %.2f\n"
                   "�����-3 � 5 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
			if(i == 4)
			printf("�����-6 � 1 ���������� = %.2f\n"
                   "�����-6 � 2 ���������� = %.2f\n"
                   "�����-6 � 3 ���������� = %.2f\n"
                   "�����-6 � 4 ���������� = %.2f\n"
                   "�����-6 � 5 ���������� = %.2f\n",
                   	recount[0][i], recount[1][i], recount[2][i], recount[3][i], recount[4][i]);
        }
    	printf("������ ������ � 100� �������� = %.2f\n", ratio_calc[i]);
    	puts("=================================================");
    }

    printf("����������� ��������������� ������������ (��� i=3) = %.2f\n", result1);
    printf("����������� ��������������� ������������ (��� i=5) = %.2f\n", result2);
    // �������� ��� �� ��� �� 0 �� 1
	puts("=================================================");
    printf("������ ��������� �������? ������� 0 ���� ��� � ����� ������ ����� ���� ��:\n>");
    scanf(" %d", &ans);
    fflush(stdin);
    }
    while(ans);
    return 0;
}
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/*�����*/
// ���������� sign ��� ������ �������, ����� ������ ��������� ���������
// ����������� �� ����� � �������� ������
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
	double numerator = 0.0; // ��������� � �������
	double denominator = 0.0; // ����������� � �������
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
/*������*/
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
