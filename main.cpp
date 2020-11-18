

//����1����Ȼ���׻�������������������ǲ�ȷ���Ƿ�����ǡ��װ�����������
//����2����������ͨ�����ϵĵ������ﵽ���õ�Ч��������Ҳ�п��ܳ������յ�����
//����3����������������ѡ��װ�䲢δ�������ǣ���Ҫͨ������Ч���ж����ѡ��
//����4����������ʵ���ݴ����������⣨�߶ȹ��ߵı��䣬�Լ�װ�������⣩

// #include <jni.h>
#include"Point.h"
#include"Cargo.h"
#include"Chromosome.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include "json.hpp"
#include <string>
#include <set>
#include <cstdint>
#include<ctime>

template <typename Container>
void delete_vector(const Container& container, int size)
{
    for (int i = 0; i < size; ++i)
        delete container[i];
}

template <typename Container>
void delete_vector(const Container& container)
{
    for (auto ptr : container)
        delete ptr;
}

using namespace std;
using json = nlohmann::json;
//���Ӿ�̬��
#pragma comment( lib, "ws2_32.lib" )
/*
UDP��
��ʼ�������
1.�����׽���
2.���׽���
3.��������
4.�ر�����
*/
//DWORD WINAPI RecvFun(LPVOID lp);
int SUM = 0;        //��ʾ�Ŵ��㷨��������
int CHE = 0;    //��ʾ��ʾ���ı��
int ONE = 1;    //��ʾһ��һ����ʾ
int minCargo1 = 2000; //�����ռ��С��ems(�ֱ��ʾ��һ�ڶ�����С)
int minCargo2 = 2000;
int minCargo3 = 2000;

//================================================��װ���Լ������ʼ��Ϣ
Point p1_(0, 0, 0);
Point p2_(720, 230, 215);

string global_tradeId;

Bin* bin[maxnbin_]; //��װ�������
Cargo *ca[maxn_];   //���������
//double SUM[maxnbin_];   //ÿ������װ����
//int LengthOfPacking = 100;
/*
void init()
{
//��װ��Ĳ���
    for(int i = 0; i < maxnbin; i++)
    {
        bin[i] = new Bin(p1_, p2_);
        bin[i]->SetNum(i);
    }
}
*/



//�������Ϣ�洢
vector<Cargo*>vecc;
vector<Cargo*>::iterator itc;
//�����пռ�����
vector<EMSpace>vece;
//����ά�������пռ�
vector<EMSpace>vece2;
vector<EMSpace>::iterator ite;
vector<EMSpace>::iterator ite2;

//================================================�Ŵ��㷨�����Լ�����
const int MAXN = 1;    //��ʾ�Ŵ��㷨ÿһ����Ⱥ������
const int EV = 0; //��ʾ�Ŵ��㷨�����Ĵ���
const double SaProbability = 0.100; //����Ⱦɫ�屣���ĸ���
const double OvProbability = 0.600; //����ĸ���
const double VaProbability = 0.300;  //����ĸ���
const double OvElite = 0.500;   //�ӽ�����ѡ��Ӣ����
const double VaGene = 0.400;    //����Ⱦɫ�����ĸ���

vector<Chromosome*>vec_now;    //������Ⱥ
vector<Chromosome*>vec_next;    //������һ����Ⱥ
vector<Chromosome*>vec_elite;   //��ӢȾɫ����Ⱥ
Chromosome* chrom[MAXN];        //Ⱦɫ�����飨��������ɣ�

vector<Chromosome*> allocated_chrosomes;

void Init()
{
    vec_now.clear();
    vec_next.clear();
    vec_elite.clear();
    vecc.clear();
    vece.clear();
    vece2.clear();
    allocated_chrosomes.clear();
}

void Initialize();  //�����ʼ����Ⱥ�� �õ���һ������
void CaculaFitness();   //����������Ӧֵ

void Selecting();     //��Ⱥѡ��
void Crossing();    //�ӽ�
void Variating();   //����
Chromosome* GeneticAlgorithm();    //�Ŵ��㷨(��󷵻�һ����ѵ�Ⱦɫ��)

//==============================================װ���㷨�����Լ�����

void ShowPoint(Point p1, Point p2);     //��ʾ�������������
int Compare(Point p1, Point p2);        //�Ƚ��������λ�ã� ����ֵ1��ʾp1��p2���Һ��Ϸ�
bool  DFTRC_2(vector<Cargo*>vec, Cargo *c, Bin *b, EMSpace*);      //���ӵı�š�װ��ռ��ѡ��
// EMSpace DFTRC_2(vector<Cargo*>vec, Cargo *c, Bin *b);      //���ӵı�š�װ��ռ��ѡ��

int CargoOrientation(vector<Cargo*>vec, Cargo *c, EMSpace ems);    //װ�䷽���ѡ��
void CargoPacking(Cargo *c, EMSpace ems, int flag);   //����װ������ӣ� װ��Ŀռ䣬 ���ӵķ�������Ԫ��
void StateUpdate(Cargo *c, int flag);         //ά��װ��Ŀռ�
void ShowCargo(vector<Cargo*> vec);     //װ������λ�����
//void SocketUDP(const char* sendBuf);       //ͨ��UDP����ͨ��
//����UDPͨ�Ź���
void SocketUDP(const char* sendBuf)
{
    // WSADATA wsaData;
    // //��ʼ��socket
    // WSAStartup(MAKEWORD(2,2),&wsaData);
    // //�׽��ֳ�ʼ��
    // SOCKET sendSocket;
    // sendSocket=socket(AF_INET,SOCK_DGRAM,IPPROTO_UDP);
    // //���ü�������������ַ
    // sockaddr_in seAddr;
    // seAddr.sin_family=AF_INET;
    // seAddr.sin_port=htons(55554);
    // seAddr.sin_addr.S_un.S_addr=inet_addr("127.0.0.1");//;htonl(INADDR_ANY)

    // //��ʼ��
    // int bufLen=strlen(sendBuf);
    // //��
    // sendto(sendSocket,sendBuf,bufLen,0,(SOCKADDR *)&seAddr,sizeof(seAddr));
    // //������ɣ��ر�socket
    // closesocket(sendSocket);
    // //�ͷ���Դ���˳�
    // WSACleanup();
}
//===============================================�Ŵ��㷨��ķ�����ʵ��
bool ComparePacking1(Cargo* x, Cargo* y)     //�����Ӱ���Ⱦɫ������������
{
    if(x->GetBPS() != y->GetBPS())
        return x->GetBPS() < y->GetBPS();
    else
        return x->GetId() < y->GetId();
}
bool ComparePacking2(Cargo* x, Cargo* y)    //�����Ӱ�����Ž�������
{
    return x->GetId() < y->GetId();
}
double Chromosome::Fitness(vector<Cargo*> vec, Bin* b)        //Ⱦɫ�����Ӧֵ���㺯��
{
    //��װ���ʼ��Ϣά��
    vece.clear();
    EMSpace ems(b->GetPoint1(), b->GetPoint2());
    vece.push_back(ems);
    sort(vec.begin(), vec.end(), ComparePacking2);

    int k = 0;  //��ʾȾɫ�����ֵ����
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        (*it)->SetBPS(ch[k] + (*it)->ding);
        (*it)->SetVBO(ch[maxn + k]);
        k++;
    }

    sort(vec.begin(), vec.end(), ComparePacking1);

    int numk = 0;   //��ʾ��װ��ı��
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        //// << (*it)->GetLength() << " " << (*it)->GetWidth() << " " << (*it)->GetHigh() << endl;
        // EMSpace ems_select = DFTRC_2(vec, *it, bin[numk]);   //ѡ������װ��Ŀռ�
        EMSpace ems_select(ems);
        bool valid = DFTRC_2(vec, *it, bin[numk], &ems_select);   //ѡ������װ��Ŀռ�
        if (!valid)
            return 0;

        //// << ems_select.GetLength() << " " << ems_select.GetWidth() << " " << ems_select.GetHigh() << endl;
        numk = (*it)->numofpack;
        //// << numk;
        int flag = CargoOrientation(vec, *it, ems_select); //���ؼ�Ϊ���ӵķ���
        //// << flag << endl;
        CargoPacking(*it, ems_select, flag);  //����װ��
        StateUpdate(*it, flag);
    }

    {
        // keep suite?

        map<string, map<int, map<string, int>>> table;

        for (Cargo* cargo : vec)
        {
            const string& set_code = cargo->set_code;
            if (set_code != "")
            {
                table[set_code][cargo->numofpack][cargo->materialCode]++;
            }
        }

        for (auto& same_set : table)
        {
            for (auto&& same_bin : same_set.second)
            {
                const map<string, int>& cargo_amounts = same_bin.second;
                int amount = cargo_amounts.begin()->second;

                for (auto&& pair : cargo_amounts)
                {
                    if (pair.second != amount)
                    {
                        fitness = 0;
                        return 0;
                    }
                }
            }
        }
    }

    vece.clear();//װ��ռ�������
    int max_x = 0, max_y = 0, max_z = 0;

    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        if((*it)->numofpack == numk)
        {
            max_x = max(max_x, ((*it)->GetPoint2()).GetPointX());
            max_y = max(max_y, ((*it)->GetPoint2()).GetPointY());
            max_z = max(max_z, ((*it)->GetPoint2()).GetPointZ());
        }
    }
    double value_k = numk + (((max_x) * 1.0) / ((bin[numk]->GetPoint2()).GetPointX())* 1.0);  //��һ������������ʾװ��һ����Ҫ���ٸ�����
    //double value_k = numk + (((max_x * max_y * max_z) * 1.0) / (p2_.GetPointX() * p2_.GetPointY() * p2_.GetPointZ())* 1.0);  //��һ������������ʾװ��һ����Ҫ���ٸ�����
    value_k = 1.0 / value_k;
    fitness = value_k;
    sort(vec.begin(), vec.end(), ComparePacking2);


    return value_k;
}

bool Chromosome::PackingOf(vector<Cargo*> vec)
{

    // printf("װ����ʾ\n");
    //��װ���ʼ��Ϣά��
    vece.clear();
    EMSpace ems(bin[0]->GetPoint1(), bin[0]->GetPoint2());
    vece.push_back(ems);
    sort(vec.begin(), vec.end(), ComparePacking2);
    int k = 0;  //��ʾȾɫ�����ֵ����
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        int ding = (*it)->ding;
        (*it)->SetBPS(ch[k] + ding);
        (*it)->SetVBO(ch[maxn + k]);
        k++;
    }
    sort(vec.begin(), vec.end(), ComparePacking1);

    int numk = 0;   //��ʾ��װ��ı��                                          //��ʾװ�����ĸ���װ��
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        // EMSpace ems_select = DFTRC_2(vec, *it, bin[numk]);   //ѡ������װ��Ŀռ�
        EMSpace ems_select(ems);
        bool valid = DFTRC_2(vec, *it, bin[numk], &ems_select);   //ѡ������װ��Ŀռ�
        if (!valid)
            return false;

        numk = (*it)->numofpack;
        // // << numk;
        int flag = CargoOrientation(vec, *it, ems_select); //���ؼ�Ϊ���ӵķ���

        (*it)->flag = flag;
        // if((*it)->flag == 2)
        //     (*it)->flag = 3;
        // else if((*it)->flag == 3)
        //     (*it)->flag = 2;
        //// << flag << endl;
        CargoPacking(*it, ems_select, flag);  //����װ��
        StateUpdate(*it, flag);
        //sum_all += (*it)->GetLength() * (*it)->GetWidth() * (*it)->GetHigh();
    }
    vece.clear();//װ��ռ�������
    //��װ���˳���λ�÷���

    sort(vec.begin(), vec.end(), ComparePacking1);

    double SUM[maxnbin_];
    memset(SUM, 0, sizeof(SUM));
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        int np = (*it)->numofpack;
        SUM[np] += (*it)->GetLength() * (*it)->GetWidth() * (*it)->GetHigh();
    }
    for(int i=0; i<=numk; i++)
    {
        cerr << "packing rate of car " << i << " = " << SUM[i] * 1.0 / (bin[i]->GetLength() * bin[i]->GetWidth() * bin[i]->GetHigh()) * 1.0 << endl;
    }
    // ShowCargo(vec);

    return true;
}
//=======================================================�Ŵ��㷨������ʵ��

void Initialize()   //�����ʼ����Ⱥ�� �õ���һ����Ⱥ
{
    vec_now.clear();
    srand(time(NULL));  //�������������
    ID = 0;
    for(int i=0; i<MAXN; i++)
    {
        double ch_variable[maxn * 2];
        for(int j=0; j<maxn * 2; j++)
        {
            double p = (rand() % 1000)/ 1000.0;
            ch_variable[j] = p;
        }/*
        // for(int j = 0; j<maxn * 2; j++)
        // {
        //     // << ch_variable[j] << " ";
        // }
        // // << endl;
        */
        chrom[i] = new Chromosome(ch_variable);
        allocated_chrosomes.push_back(chrom[i]);
        chrom[i]->SetId(ID);
        ID++;
        vec_now.push_back(chrom[i]);
    }
    // // << "Ⱦɫ��������С:" << vec_now.size() << endl;
}
void CaculaFitness()    //������Ӧֵ����ÿһ������Ⱦɫ�嶼������
{
    double sum = 0; //��Ӧֵ�����ܺ�
    double value[vec_now.size()];
    int k = 0;  //���ڶ�valueֵ���б���
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)
    {
        //if((*it)->GetFitness() == -1)
        value[k] = (*it)->Fitness(vecc, bin[0]);    //�ӵ�0�����ӿ�ʼ
        //else
        //  value[k] = (*it)->GetFitness();
        sum += value[k];
        k++;
    }
    k = 0;
    double max_ = 0;
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)
    {
        double m = value[k] * 1.0/ sum * 1.0;
        (*it)->SetReFitness(m);     //����Խ��Խ��
        //// << 1.0 / (*it)->GetFitness() * 1.0 << endl;
        max_ = max((*it)->GetFitness(), max_);                                                                 //���ڵ�����ʾ��Ҫ�����㹻��Ĵ���
        k++;
    }
    SUM++;
    // // << SUM << " " << max_ << endl;
}
bool CompareChromosome(Chromosome* x, Chromosome* y)    //Ⱦɫ�������㷨
{
    if(x->GetReFitness() != y->GetReFitness())
        return x->GetReFitness() > y->GetReFitness();
    else
        return x->GetId() < y->GetId();
}
//ѡ��ӢȾɫ�屣��
void Selecting()
{
    vec_next.clear();
    vec_elite.clear();

    sort(vec_now.begin(), vec_now.end(), CompareChromosome);
    int Length = (int)(vec_now.size() * SaProbability);    //��ӢȾɫ��ĸ���
    int k = 0;
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it != vec_now.end(); it++)
    {
        if(k < Length)
        {
            vec_next.push_back(*it);        //����ӢȾɫ��ѡ����һ����
            vec_elite.push_back(*it);
            //// << "��Ӣ" << (*it)->GetReFitness() << endl << endl;
        }
        else
            break;
        k++;
    }
    //���̶��㷨 ����Ⱦɫ�����ֵ
    double sum = 0;
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)
    {
        sum = sum + (*it)->GetReFitness();
        (*it)->SetSumFitness(sum);
        //// <<(*it)->GetSumFitness() << endl;
    }
    //// << vec_next.size() << endl;
}
//�ӽ�����Ⱦɫ��
void Crossing()
{
    //// << "Test2" << endl;
    srand(time(NULL));
    int Length = vec_elite.size();    //��ӢȾɫ��ĸ���
    for(int j = 0; j < (int)(vec_now.size() * OvProbability); j++)  //�ӽ��Ĺ���
    {
        Chromosome* chrom1 = vec_now.front(), *chrom2 = vec_now.front();
        int k_ans = rand()%Length; //��ѡ��ľ�ӢȾɫ��ı��
        int k_now = 0;  //��ʾ���ڱ���Ⱦɫ��ı��
        for(vector<Chromosome*>::iterator it = vec_elite.begin(); it!= vec_elite.end(); it++)//ѡ���һ��Ⱦɫ��
        {
            if(k_now == k_ans)
            {
                chrom1 = *it;
                break;
            }
            k_now++;
        }
        double value = (rand()%1000 + 1)/1000.0;    //���̶�ѡ��ĸ���ֵ
        for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)       //ѡ��ڶ���Ⱦɫ��
        {
            if(value <= (*it)->GetSumFitness() && value > ((*it)->GetSumFitness() - (*it)->GetReFitness()))
            {
                chrom2 = *it;
                break;
            }
        }
        if(chrom1 == chrom2)
        {
            j--;
        }
        else
        {
            //// << chrom1->GetReFitness() << " " << chrom2->GetReFitness() << endl;
            double ans[maxn * 2];
            for(int i=0; i<maxn * 2; i++)
            {
                value = (rand()%1000 + 1)/1000.0;    //���̶�ѡ��ĸ���ֵ
                if(value > 0 && value <=OvElite)
                {
                    ans[i] = chrom1->GetValue(i);
                }
                else if(value > OvElite)
                {
                    ans[i] = chrom2->GetValue(i);
                }
            }
            Chromosome* chrom_ans = new Chromosome(ans);
            allocated_chrosomes.push_back(chrom_ans);
            vec_next.push_back(chrom_ans);
        }
    }
    //// << vec_next.size() << endl;
}
//�������Ⱦɫ��
void Variating()
{
    //// << "Test3" << endl;
    srand(time(NULL));
    int Length = vec_now.size() - vec_next.size();  //����Ⱦɫ�������
    //// << Length << endl;
    for(int i=0; i<Length; i++)
    {
        Chromosome* chrom = vec_now.front();   //ѡ������Ⱦɫ��
        double value = (rand()%1000 + 1)/1000.0;    //���̶�ѡ��ĸ���ֵ
        for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)       //ѡ�����Ⱦɫ��
        {
            if(value <= (*it)->GetSumFitness() && value > ((*it)->GetSumFitness() - (*it)->GetReFitness()))
            {
                chrom = *it;
                break;
            }
        }
        double ans[maxn * 2];
        for(int i=0; i<maxn * 2; i++)
        {
            value = (rand()%1000 + 1)/1000.0;    //���̶�ѡ��ĸ���ֵ
            if(value > 0 && value <=VaGene)
            {
                ans[i] = (rand()%1000 + 1)/1000.0;
            }
            else if(value > VaGene)
            {
                ans[i] = chrom->GetValue(i);
            }
        }
        Chromosome* chrom_ans = new Chromosome(ans);
        allocated_chrosomes.push_back(chrom_ans);
        vec_next.push_back(chrom_ans);
    }
    vec_now.clear();
    vec_elite.clear();
    ID = 0;
    for(vector<Chromosome*>::iterator it = vec_next.begin(); it != vec_next.end(); it++)
    {
        (*it)->SetId(ID);
        ID++;
        vec_now.push_back(*it);
    }
    vec_next.clear();
    //// << vec_now.size() << endl;
}

Chromosome* GeneticAlgorithm() //�Ŵ��㷨
{
    Initialize();   //��ʼ����Ⱥ��������ɵ�һ������

    for(int i = 0; i < EV; i++) //ÿһ�������Ŵ��㷨
    {
        CaculaFitness();    //����ÿ��Ⱦɫ����Ӧֵ����
        Selecting();    //ѡ��Ⱦɫ��
        Crossing();     //�ӽ�Ⱦɫ��
        Variating();    //����Ⱦɫ��
    }

    CaculaFitness();   //�������һ��Ⱦɫ�����Ӧֵ����

    double max_ = 0;
    Chromosome* ch = nullptr; //���ձ�ѡ���Ⱦɫ��
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it != vec_now.end(); it++)
    {
        //// << (*it)->GetFitness() << endl;
        if((*it)->GetReFitness() > max_)
        {
            max_ = (*it)->GetReFitness();
            ch = *it;
        }
    }
    return ch;
}
//==========================================================װ�����
//������������ʾ
void ShowPoint(Point p1, Point p2)
{
    // // << p1.GetPointX() << " " << p1.GetPointY() << " " << p1.GetPointZ() << endl;
    // // << p2.GetPointX() << " " << p2.GetPointY() << " " << p2.GetPointZ() << endl;
}
int Compare(Point p1, Point p2)
{
    int k = 0;
    if(p1.GetPointX() >= p2.GetPointX() && p1.GetPointY() >= p2.GetPointY() && p1.GetPointZ() >= p2.GetPointZ())
    {
        k = 1;   //˵����1�ڵ�2���Һ��ϲ�
    }
    if (p1.GetPointX() <= p2.GetPointX() && p1.GetPointY() <= p2.GetPointY() && p1.GetPointZ() <= p2.GetPointZ())
    {
        k = -1;   //˵����2�ڵ�1���Һ��ϲ�
    }
    if(p1.GetPointX() == p2.GetPointX() && p1.GetPointY() == p2.GetPointY() && p1.GetPointZ() == p2.GetPointZ())
    {
        k = 2;   //˵����1�͵�2��ͬ
    }
    return k;
}

//ѡ������װ��Ŀռ�(���������Դ�������)
// EMSpace DFTRC_2(vector<Cargo*>vec, Cargo *c, Bin *b)
bool DFTRC_2(vector<Cargo*>vec, Cargo *c, Bin *b, EMSpace* space)
{
    //// << endl;
    int k = b->GetNum();    // ���ӱ��
    double maxdist = -1;
    vector<EMSpace>::iterator it_ems = vece.end();
    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        //// <<endl;
        int x = (ite->GetPoint1()).GetPointX();
        int y = (ite->GetPoint1()).GetPointY();
        int z = (ite->GetPoint1()).GetPointZ();
        //����1
        //if(ite->GetLength() >= c->GetLength() && ite->GetWidth() >= c->GetWidth() && ite->GetHigh() >= c->GetHigh())
        if((c->node[0]).f == 1 &&c->weight <= ite->weight && ite->GetLength() >= c->GetLength() && ite->GetWidth() >= c->GetWidth() && ite->GetHigh() >= c->GetHigh() && (c->node[0]).l >= ite->level && ite->yes == 1)
        {
            if (((ite->cargo_code == c->materialCode) && (ite->layer > 0)) || ((ite->cargo_code != c->materialCode) && (ite->max_weight > c->weight)))
            {
                int k = 0;  //��k = 5��ʾ���Էŵ���
                Point point1(x, y, z);
                Point point2(x + c->GetLength(), y + c->GetWidth(), z + c->GetHigh());  //װ�������
                if(z == 0)
                {
                    k = 5;
                }
                else    //�ѻ�����
                {
                    Point *point[5];
                    point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                    point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                    {
                        if(*it == c)
                        {
                            break;
                        }
                        Point point1_ = (*it)->GetPoint1();
                        Point point2_ = (*it)->GetPoint2();
                        for(int i = 0; i < 5; i++)
                        {
                            if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                            {
                                if(a[i] == 0)
                                {
                                    a[i] = 1;
                                    k++;
                                }
                            }
                        }
                        if(k == 5)
                            break;
                    }

                    delete_vector(point);
                }

                double dist =pow((b->GetLength() - x - c->GetLength())*1.0, 2.2)  + pow((b->GetWidth() - y - c->GetWidth())*1.0, 2.0) + pow((b->GetHigh() - z - c->GetHigh())*1.0, 2.0);

                if(dist > maxdist && k == 5)
                {
                    maxdist = dist;
                    it_ems = ite;
                    //// << "        1" << endl;
                }
            }
        }

        //����2
        //if(ite->GetLength() >= c->GetLength() && ite->GetWidth() >= c->GetHigh() && ite->GetHigh() >= c->GetWidth())
        if((c->node[1]).f == 1 &&c->weight <= ite->weight&& ite->GetLength() >= c->GetLength() && ite->GetWidth() >= c->GetHigh() && ite->GetHigh() >= c->GetWidth() && (c->node[1]).l >= ite->level && ite->yes == 1)
        {
            if (((ite->cargo_code == c->materialCode) && (ite->layer > 0)) || (ite->cargo_code != c->materialCode) && (ite->max_weight > c->weight))
            {
                int k = 0;  //��k = 5��ʾ���Էŵ���
                Point point1(x, y, z);
                Point point2(x + c->GetLength(), y + c->GetHigh(), z + c->GetWidth());  //װ�������
                if(z == 0)
                {
                    k = 5;
                }
                else    //�ѻ�����
                {
                    Point *point[5];
                    point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                    //// << point[0]->GetPointX() << " " << point[0]->GetPointY() << " " << point[0]->GetPointZ() << endl;
                    point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                    {
                        if(*it == c)
                        {
                            break;
                        }
                        Point point1_ = (*it)->GetPoint1();
                        Point point2_ = (*it)->GetPoint2();
                        for(int i = 0; i < 5; i++)
                        {
                            if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                            {
                                if(a[i] == 0)
                                {
                                    a[i] = 1;
                                    k++;
                                }
                            }
                        }
                        if(k == 5)
                            break;
                    }

                    delete_vector(point);
                }
                double dist =pow((b->GetLength() - x - c->GetLength())*1.0, 2.2)  + pow((b->GetWidth() - y - c->GetHigh())*1.0, 2.0) + pow((b->GetHigh() - z - c->GetWidth())*1.0, 2.0);
                if(dist > maxdist && k == 5)
                {
                    maxdist = dist;
                    it_ems = ite;
                    //// << "        2" << endl;
                }
            }
        }

        //����3
        //if(ite->GetLength() >= c->GetWidth() && ite->GetWidth() >= c->GetLength() && ite->GetHigh() >= c->GetHigh())
        if((c->node[2]).f == 1  &&c->weight <= ite->weight && ite->GetLength() >= c->GetWidth() && ite->GetWidth() >= c->GetLength() && ite->GetHigh() >= c->GetHigh() && (c->node[2]).l >= ite->level && ite->yes == 1)
        {
            if (((ite->cargo_code == c->materialCode) && (ite->layer > 0)) || (ite->cargo_code != c->materialCode) && (ite->max_weight > c->weight))
            {
                int k = 0;  //��k = 5��ʾ���Էŵ���
                Point point1(x, y, z);
                Point point2(x + c->GetWidth(), y + c->GetLength(), z + c->GetHigh());  //װ�������
                if(z == 0)
                {
                    k = 5;
                }
                else    //�ѻ�����
                {
                    Point *point[5];
                    point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                    point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                    {
                        if(*it == c)
                        {
                            break;
                        }
                        Point point1_ = (*it)->GetPoint1();
                        Point point2_ = (*it)->GetPoint2();
                        for(int i = 0; i < 5; i++)
                        {
                            if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                            {
                                if(a[i] == 0)
                                {
                                    a[i] = 1;
                                    k++;
                                }
                            }
                        }
                        if(k == 5)
                            break;
                    }

                    delete_vector(point);
                }
                double dist =pow((b->GetLength() - x - c->GetWidth())*1.0, 2.2)  + pow((b->GetWidth() - y - c->GetLength())*1.0, 2.0) + pow((b->GetHigh() - z - c->GetHigh())*1.0, 2.0);
                if(dist > maxdist && k == 5)
                {
                    maxdist = dist;
                    it_ems = ite;
                    //// << "        3" << endl;
                }
            }
        }
        //����4
        //if(ite->GetLength() >= c->GetWidth() && ite->GetWidth() >= c->GetHigh() && ite->GetHigh() >= c->GetLength())
        if((c->node[3]).f == 1 &&c->weight <= ite->weight && ite->GetLength() >= c->GetWidth() && ite->GetWidth() >= c->GetHigh() && ite->GetHigh() >= c->GetLength() && (c->node[3]).l >= ite->level && ite->yes == 1)
        {
            if (((ite->cargo_code == c->materialCode) && (ite->layer > 0)) || (ite->cargo_code != c->materialCode) && (ite->max_weight > c->weight))
            {
                int k = 0;  //��k = 5��ʾ���Էŵ���
                Point point1(x, y, z);
                Point point2(x + c->GetWidth(), y + c->GetHigh(), z + c->GetLength());  //װ�������
                if(z == 0)
                {
                    k = 5;
                }
                else    //�ѻ�����
                {
                    Point *point[5];
                    point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                    point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                    {
                        if(*it == c)
                        {
                            break;
                        }
                        Point point1_ = (*it)->GetPoint1();
                        Point point2_ = (*it)->GetPoint2();
                        for(int i = 0; i < 5; i++)
                        {
                            if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                            {
                                if(a[i] == 0)
                                {
                                    a[i] = 1;
                                    k++;
                                }
                            }
                        }
                        if(k == 5)
                            break;
                    }

                    delete_vector(point);
                }
                double dist =pow((b->GetLength() - x - c->GetWidth())*1.0, 2.2)  + pow((b->GetWidth() - y - c->GetHigh())*1.0, 2.0) + pow((b->GetHigh() - z - c->GetLength())*1.0, 2.0);
                if(dist > maxdist && k == 5)
                {
                    maxdist = dist;
                    it_ems = ite;
                    //// << "        4" << endl;
                }
            }
        }
        //����5
        //if(ite->GetLength() >= c->GetHigh() && ite->GetWidth() >= c->GetLength() && ite->GetHigh() >= c->GetWidth())
        if((c->node[4]).f == 1 &&c->weight <= ite->weight && ite->GetLength() >= c->GetHigh() && ite->GetWidth() >= c->GetLength() && ite->GetHigh() >= c->GetWidth() && (c->node[4]).l >= ite->level && ite->yes == 1)
        {
            if (((ite->cargo_code == c->materialCode) && (ite->layer > 0)) || (ite->cargo_code != c->materialCode) && (ite->max_weight > c->weight))
            {
                int k = 0;  //��k = 5��ʾ���Էŵ���
                Point point1(x, y, z);
                Point point2(x + c->GetHigh(), y + c->GetLength(), z + c->GetWidth());  //װ�������
                if(z == 0)
                {
                    k = 5;
                }
                else    //�ѻ�����
                {
                    Point *point[5];
                    point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                    point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                    {
                        if(*it == c)
                        {
                            break;
                        }
                        Point point1_ = (*it)->GetPoint1();
                        Point point2_ = (*it)->GetPoint2();
                        for(int i = 0; i < 5; i++)
                        {
                            if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                            {
                                if(a[i] == 0)
                                {
                                    a[i] = 1;
                                    k++;
                                }
                            }
                        }
                        if(k == 5)
                            break;
                    }

                    delete_vector(point);
                }
                double dist =pow((b->GetLength() - x - c->GetHigh())*1.0, 2.2)  + pow((b->GetWidth() - y - c->GetLength())*1.0, 2.0) + pow((b->GetHigh() - z - c->GetWidth())*1.0, 2.0);
                if(dist > maxdist && k == 5)
                {
                    maxdist = dist;
                    it_ems = ite;
                }
            }
        }
        //����6
        //if(ite->GetLength() >= c->GetHigh() && ite->GetWidth() >= c->GetWidth() && ite->GetHigh() >= c->GetLength())
        if((c->node[5]).f == 1 &&c->weight <= ite->weight && ite->GetLength() >= c->GetHigh() && ite->GetWidth() >= c->GetWidth() && ite->GetHigh() >= c->GetLength() && (c->node[5]).l >= ite->level && ite->yes == 1)
        {
            if (((ite->cargo_code == c->materialCode) && (ite->layer > 0)) || (ite->cargo_code != c->materialCode) && (ite->max_weight > c->weight))
            {
                int k = 0;  //��k = 5��ʾ���Էŵ���
                Point point1(x, y, z);
                Point point2(x + c->GetHigh(), y + c->GetWidth(), z + c->GetLength());  //װ�������
                if(z == 0)
                {
                    k = 5;
                }
                else    //�ѻ�����
                {
                    Point *point[5];
                    point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                    point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                    point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                    int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                    {
                        if(*it == c)
                        {
                            break;
                        }
                        Point point1_ = (*it)->GetPoint1();
                        Point point2_ = (*it)->GetPoint2();
                        for(int i = 0; i < 5; i++)
                        {
                            if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                            {
                                if(a[i] == 0)
                                {
                                    a[i] = 1;
                                    k++;
                                }
                            }
                        }
                        if(k == 5)
                            break;
                    }

                    delete_vector(point);
                }
                double dist =pow((b->GetLength() - x - c->GetHigh())*1.0, 2.2)  + pow((b->GetWidth() - y - c->GetWidth())*1.0, 2.0) + pow((b->GetHigh() - z - c->GetLength())*1.0, 2.0);
                if(dist > maxdist && k == 5)
                {
                    maxdist = dist;
                    it_ems = ite;
                }
            }
        }
    }

    if(it_ems == vece.end())
    {
        k++;

        if(!bin[k])
            return false;

        vece.clear();
        Point P1 = bin[k]->GetPoint1();
        Point P2 = bin[k]->GetPoint2();
        EMSpace ems_(P1, P2);
        vece.push_back(ems_);
        c->numofpack = bin[k]->GetNum();
        *space = ems_;
        // return ems_;
        return true;
    }
    else
    {
        c->numofpack = b->GetNum();
        // return *it_ems; //��������װ��Ŀռ�
        *space = *it_ems;
        return true;
    }
}

//ȷ��װ��ķ���
int CargoOrientation(vector<Cargo*>vec, Cargo *c, EMSpace ems)
{
    c->Setnum(0);
    for(int j=0; j<6; j++)
    {
        c->SetOri(j, 0);
    }
    //// << (c->node[0]).f << " " << c->weight << " " << ems.weight << " " << (c->node[0]).l << " " << ems.level << " " << ems.yes << endl;
    //����1
    //if(ems.GetLength() >= c->GetLength() && ems.GetWidth() >= c->GetWidth() && ems.GetHigh() >= c->GetHigh())
    if((c->node[0]).f == 1&& c->weight <= ems.weight && (c->node[0]).l >= ems.level && ems.yes == 1 && ems.GetLength() >= c->GetLength() && ems.GetWidth() >= c->GetWidth() && ems.GetHigh() >= c->GetHigh())
    {
        if (((ems.cargo_code == c->materialCode) && (ems.layer > 0)) || (ems.cargo_code != c->materialCode) && (ems.max_weight > c->weight))
        {
            //// << "yes" << endl;
            int k = 0;  //��k = 5��ʾ���Էŵ���
            int x = (ems.GetPoint1()).GetPointX();
            int y = (ems.GetPoint1()).GetPointY();
            int z = (ems.GetPoint1()).GetPointZ();

            Point point1(x, y, z);
            Point point2(x + c->GetLength(), y + c->GetWidth(), z + c->GetHigh());  //װ�������
            if(z == 0)
            {
                k = 5;
            }
            else    //�ѻ�����
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                {
                    if(*it == c)
                    {
                        break;
                    }
                    Point point1_ = (*it)->GetPoint1();
                    Point point2_ = (*it)->GetPoint2();
                    for(int i = 0; i < 5; i++)
                    {
                        if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                        {
                            if(a[i] == 0)
                            {
                                a[i] = 1;
                                k++;
                            }
                        }
                    }
                    if(k == 5)
                        break;
                }

            }
            if(k == 5)
            {
                int k_ = c->Getnum();   //�ɷ����Ӹ���ѡ��
                c->SetOri(k_, 1);
                c->Setnum(k_ + 1);
            }
        }
    }

    //����2
    //if(ems.GetLength() >= c->GetLength() && ems.GetWidth() >= c->GetHigh() && ems.GetHigh() >= c->GetWidth())
    if((c->node[1]).f == 1 && c->weight <= ems.weight &&(c->node[1]).l >= ems.level && ems.yes == 1&& ems.GetLength() >= c->GetLength() && ems.GetWidth() >= c->GetHigh() && ems.GetHigh() >= c->GetWidth())
    {
        if (((ems.cargo_code == c->materialCode) && (ems.layer > 0)) || (ems.cargo_code != c->materialCode) && (ems.max_weight > c->weight))
        {
            int k = 0;  //��k = 5��ʾ���Էŵ���
            int x = (ems.GetPoint1()).GetPointX();
            int y = (ems.GetPoint1()).GetPointY();
            int z = (ems.GetPoint1()).GetPointZ();

            Point point1(x, y, z);
            Point point2(x + c->GetLength(), y + c->GetHigh(), z + c->GetWidth());  //װ�������
            if(z == 0)
            {
                k = 5;
            }
            else    //�ѻ�����
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                {
                    if(*it == c)
                    {
                        break;
                    }
                    Point point1_ = (*it)->GetPoint1();
                    Point point2_ = (*it)->GetPoint2();
                    for(int i = 0; i < 5; i++)
                    {
                        if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                        {
                            if(a[i] == 0)
                            {
                                a[i] = 1;
                                k++;
                            }
                        }
                    }
                    if(k == 5)
                        break;
                }
            }
            if(k == 5)
            {
                int k_ = c->Getnum();   //�ɷ����Ӹ���ѡ��
                c->SetOri(k_, 2);
                c->Setnum(k_ + 1);
            }
        }
    }

    //����3
    //if(ems.GetLength() >= c->GetWidth() && ems.GetWidth() >= c->GetLength() && ems.GetHigh() >= c->GetHigh())
    if((c->node[2]).f == 1 && c->weight <= ems.weight &&(c->node[2]).l >= ems.level && ems.yes == 1&& ems.GetLength() >= c->GetWidth() && ems.GetWidth() >= c->GetLength() && ems.GetHigh() >= c->GetHigh())
    {
        if (((ems.cargo_code == c->materialCode) && (ems.layer > 0)) || (ems.cargo_code != c->materialCode) && (ems.max_weight > c->weight))
        {
            int k = 0;  //��k = 5��ʾ���Էŵ���
            int x = (ems.GetPoint1()).GetPointX();
            int y = (ems.GetPoint1()).GetPointY();
            int z = (ems.GetPoint1()).GetPointZ();

            Point point1(x, y, z);
            Point point2(x + c->GetWidth(), y + c->GetLength(), z + c->GetHigh());  //װ�������
            if(z == 0)
            {
                k = 5;
            }
            else    //�ѻ�����
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                {
                    if(*it == c)
                    {
                        break;
                    }
                    Point point1_ = (*it)->GetPoint1();
                    Point point2_ = (*it)->GetPoint2();
                    for(int i = 0; i < 5; i++)
                    {
                        if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                        {
                            if(a[i] == 0)
                            {
                                a[i] = 1;
                                k++;
                            }
                        }
                    }
                    if(k == 5)
                        break;
                }
            }
            if(k == 5)
            {
                int k_ = c->Getnum();   //�ɷ����Ӹ���ѡ��
                c->SetOri(k_, 3);
                c->Setnum(k_ + 1);
            }
        }
    }
    //����4
    //if(ems.GetLength() >= c->GetWidth() && ems.GetWidth() >= c->GetHigh() && ems.GetHigh() >= c->GetLength())
    if((c->node[3]).f == 1 && c->weight <= ems.weight &&(c->node[3]).l >= ems.level && ems.yes == 1 && ems.GetLength() >= c->GetWidth() && ems.GetWidth() >= c->GetHigh() && ems.GetHigh() >= c->GetLength())
    {
        if (((ems.cargo_code == c->materialCode) && (ems.layer > 0)) || (ems.cargo_code != c->materialCode) && (ems.max_weight > c->weight))
        {
            int k = 0;  //��k = 5��ʾ���Էŵ���
            int x = (ems.GetPoint1()).GetPointX();
            int y = (ems.GetPoint1()).GetPointY();
            int z = (ems.GetPoint1()).GetPointZ();

            Point point1(x, y, z);
            Point point2(x + c->GetWidth(), y + c->GetHigh(), z + c->GetLength());  //װ�������
            if(z == 0)
            {
                k = 5;
            }
            else    //�ѻ�����
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                {
                    if(*it == c)
                    {
                        break;
                    }
                    Point point1_ = (*it)->GetPoint1();
                    Point point2_ = (*it)->GetPoint2();
                    for(int i = 0; i < 5; i++)
                    {
                        if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                        {
                            if(a[i] == 0)
                            {
                                a[i] = 1;
                                k++;
                            }
                        }
                    }
                    if(k == 5)
                        break;
                }
            }
            if(k == 5)
            {
                int k_ = c->Getnum();   //�ɷ����Ӹ���ѡ��
                c->SetOri(k_, 4);
                c->Setnum(k_ + 1);
            }
        }
    }
    //����5
    //if(ems.GetLength() >= c->GetHigh() && ems.GetWidth() >= c->GetLength() && ems.GetHigh() >= c->GetWidth())
    if((c->node[4]).f == 1 && c->weight <= ems.weight &&(c->node[4]).l >= ems.level && ems.yes == 1 && ems.GetLength() >= c->GetHigh() && ems.GetWidth() >= c->GetLength() && ems.GetHigh() >= c->GetWidth())
    {
        if (((ems.cargo_code == c->materialCode) && (ems.layer > 0)) || (ems.cargo_code != c->materialCode) && (ems.max_weight > c->weight))
        {
            int k = 0;  //��k = 5��ʾ���Էŵ���
            int x = (ems.GetPoint1()).GetPointX();
            int y = (ems.GetPoint1()).GetPointY();
            int z = (ems.GetPoint1()).GetPointZ();

            Point point1(x, y, z);
            Point point2(x + c->GetHigh(), y + c->GetLength(), z + c->GetWidth());  //װ�������
            if(z == 0)
            {
                k = 5;
            }
            else    //�ѻ�����
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                {
                    if(*it == c)
                    {
                        break;
                    }
                    Point point1_ = (*it)->GetPoint1();
                    Point point2_ = (*it)->GetPoint2();
                    for(int i = 0; i < 5; i++)
                    {
                        if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                        {
                            if(a[i] == 0)
                            {
                                a[i] = 1;
                                k++;
                            }
                        }
                    }
                    if(k == 5)
                        break;
                }
            }
            if(k == 5)
            {
                int k_ = c->Getnum();   //�ɷ����Ӹ���ѡ��
                c->SetOri(k_, 5);
                c->Setnum(k_ + 1);
            }
        }
    }
    //����6
    //if(ems.GetLength() >= c->GetHigh() && ems.GetWidth() >= c->GetWidth() && ems.GetHigh() >= c->GetLength())
    if((c->node[5]).f == 1 && c->weight <= ems.weight &&(c->node[5]).l >= ems.level && ems.yes == 1 && ems.GetLength() >= c->GetHigh() && ems.GetWidth() >= c->GetWidth() && ems.GetHigh() >= c->GetLength())
    {
        if (((ems.cargo_code == c->materialCode) && (ems.layer > 0)) || (ems.cargo_code != c->materialCode) && (ems.max_weight > c->weight))
        {
            int k = 0;  //��k = 5��ʾ���Էŵ���
            int x = (ems.GetPoint1()).GetPointX();
            int y = (ems.GetPoint1()).GetPointY();
            int z = (ems.GetPoint1()).GetPointZ();

            Point point1(x, y, z);
            Point point2(x + c->GetHigh(), y + c->GetWidth(), z + c->GetLength());  //װ�������
            if(z == 0)
            {
                k = 5;
            }
            else    //�ѻ�����
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //����жϷ�
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //��ʾ����δ֧��

                for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
                {
                    if(*it == c)
                    {
                        break;
                    }
                    Point point1_ = (*it)->GetPoint1();
                    Point point2_ = (*it)->GetPoint2();
                    for(int i = 0; i < 5; i++)
                    {
                        if(point2_.GetPointZ() == point[i]->GetPointZ() && point[i]->GetPointX() >= point1_.GetPointX() && point[i]->GetPointX() <= point2_.GetPointX() && point[i]->GetPointY() >= point1_.GetPointY() && point[i]->GetPointY() <= point2_.GetPointY())
                        {
                            if(a[i] == 0)
                            {
                                a[i] = 1;
                                k++;
                            }
                        }
                    }
                    if(k == 5)
                        break;
                }
            }
            if(k == 5)
            {
                int k_ = c->Getnum();   //�ɷ����Ӹ���ѡ��
                c->SetOri(k_, 6);
                c->Setnum(k_ + 1);
            }
        }
    }

    int flag = (int)(c->GetVBO() * c->Getnum());
    return c->GetOri(flag);
}

//����װ�����
void CargoPacking(Cargo *c, EMSpace ems, int flag)   //����װ������ӣ� װ��Ŀռ䣬 ���ӵķ�������Ԫ��
{
    Point cargop1 = ems.GetPoint1();
    c->SetPoint1(cargop1);
    int x_ = (ems.GetPoint1()).GetPointX();
    int y_ = (ems.GetPoint1()).GetPointY();
    int z_ = (ems.GetPoint1()).GetPointZ();

    //// << "test 1" << endl;
    //// << x_ << " " << y_ << " " << z_ << endl;
    //// << "test 2" << endl;
    //// << c->GetLength() << " " << c->GetWidth() << " " << c->GetHigh() << endl;
    //// << "test 3" << endl;
    if(flag == 1)
    {
        int x = x_ + c->GetLength();
        int y = y_ + c->GetWidth();
        int z = z_ + c->GetHigh();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //// << x << " " << y << " " << z << endl;
    }
    else if(flag == 2)
    {
        int x = x_ + c->GetLength();
        int y = y_ + c->GetHigh();
        int z = z_ + c->GetWidth();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //// << x << " " << y << " " << z << endl;
    }
    else if(flag == 3)
    {
        int x = x_ + c->GetWidth();
        int y = y_ + c->GetLength();
        int z = z_ + c->GetHigh();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //// << x << " " << y << " " << z << endl;
    }
    else if(flag == 4)
    {
        int x = x_ + c->GetWidth();
        int y = y_ + c->GetHigh();
        int z = z_ + c->GetLength();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //// << x << " " << y << " " << z << endl;
    }
    else if(flag == 5)
    {
        int x = x_ + c->GetHigh();
        int y = y_ + c->GetLength();
        int z = z_ + c->GetWidth();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //// << x << " " << y << " " << z << endl;
    }
    else if(flag == 6)
    {
        int x = x_ + c->GetHigh();
        int y = y_ + c->GetWidth();
        int z = z_ + c->GetLength();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //// << x << " " << y << " " << z << endl;
    }

}

int map_input_orientation(int orientation)
{
    return orientation;
    switch (orientation)
    {
    case 1:
        return 4;
    case 2:
        return 3;
    case 3:
        return 2;
    case 4:
        return 5;
    case 5:
        return 6;
    case 6:
        return 1;
    }
    throw runtime_error(
        string("unknown input orientation number ") + to_string(orientation));
}

int map_output_orientation(int orientation)
{
    return orientation;
    switch (orientation)
    {
    case 1:
        return 6;
    case 2:
        return 3;
    case 3:
        return 2;
    case 4:
        return 1;
    case 5:
        return 4;
    case 6:
        return 5;
    }
    throw runtime_error(string("unknown output orientation number ") + to_string(orientation));
}

int temp_convert(int orientation)
{
    if (orientation == 2)
        return 3;
    else if (orientation == 3)
        return 2;
    else
        return orientation;
}

json calcGoodList(vector<Cargo*> vec, int& index)
{
    json cargoes = json::array();

    for (Cargo* cargo : vec)
    {
        Point bottom_left = cargo->GetPoint1();
        cargoes.push_back(
        {
            {"materialCode", cargo->materialCode},
            {"restrictionFlag", map_output_orientation(temp_convert(cargo->flag)) },
            {"z", bottom_left.GetPointX() * 10 },
            {"x", bottom_left.GetPointY() * 10},
            {"y", bottom_left.GetPointZ() * 10},

            // {"x", bottom_left.GetPointX() * 10 },
            // {"y", bottom_left.GetPointY() * 10},
            // {"z", bottom_left.GetPointZ() * 10},
            {"trainIndex", index++}
        });
    }

    return cargoes;
}

json calcGoodList(const Cargo* cargo, int& index)
{
    json cargoes = json::array();

    Point bottom_left = cargo->GetPoint1();
    cargoes.push_back(
    {
        {"materialCode", cargo->materialCode},
        {"restrictionFlag", map_output_orientation(temp_convert(cargo->flag)) },
        {"z", bottom_left.GetPointX() * 10 },
        {"x", bottom_left.GetPointY() * 10},
        {"y", bottom_left.GetPointZ() * 10},

        // {"x", bottom_left.GetPointX() * 10 },
        // {"y", bottom_left.GetPointY() * 10},
        // {"z", bottom_left.GetPointZ() * 10},
        {"trainIndex", index++}
    });

    return cargoes;
}

json calcStepList(const map<string, vector<Cargo*>>& orders, int& index)
{
    int step_index = 1;
    json steps = json::array();;

    // for (auto&& pair : orders) {
    // 	steps.push_back({
    // 			{"step", step_index++},
    // 			{"qty", pair.second.size()},
    // 			{"directionNum", "1*1*1"},
    // 			{"orderCode", pair.first},
    // 			{"goodList", calcGoodList(pair.second, index)}
    // 		});
    // }

    vector<pair<string, vector<Cargo*>>> orders_vec(orders.begin(), orders.end());
    sort(orders_vec.begin(), orders_vec.end(),
         [](const pair<string, vector<Cargo*>>& px, const pair<string, vector<Cargo*>>& py)
    {
        return px.second[0]->ding < py.second[0]->ding;
    });

    for (auto&& pair : orders_vec)
    {
        for (const Cargo* cargo : pair.second)
        {
            steps.push_back(
            {
                {"step", step_index++},
                {"qty", 1},
                {"directionNum", "1*1*1"},
                {"orderCode", pair.first},
                {"goodList", calcGoodList(cargo, index) }
            });
        }
    }

    return steps;
}

int calcGoodNum(const map<string, vector<Cargo*>>& orders)
{
    int num = 0;

    for (auto&& pair : orders)
    {
        num += pair.second.size();
    }

    return num;
}

json calcTrains(const map<int, map<string, vector<Cargo*>>>& containers)
{
    json trainList = json::array();

    for (auto&& pair : containers)
    {
        int index =0;

        int64_t totalCapacityCM3 = 0;
        double totalWeight = 0;

        for (auto&& p2 : pair.second)
        {
            for (Cargo* cargo : p2.second)
            {
                totalCapacityCM3 += cargo->GetVolumeCM3();
                totalWeight += cargo->GetWeight();
            }
        }

        double totalCapacityM3 = (double)totalCapacityCM3 / 10e6;

        Bin* currentBin = bin[pair.first];

        double packingRate = (double)totalCapacityCM3 / currentBin->GetVolumeCM3();

        trainList.push_back(
        {
            { "train", pair.first},
            {"modelCode", bin[pair.first]->modelCode },
            {"goodNum", calcGoodNum(pair.second)},
            {"totalCapacity", totalCapacityM3 },  // !!!
            {"totalWeight", totalWeight}, /// !!!
            {"packingRate", packingRate}, /// !!!s
            {"stepList", calcStepList(pair.second, index) }
        });
    }

    return trainList;
}

string another_ANS(vector<Cargo*>  vec, bool valid)
{
    sort(vec.begin(), vec.end(), ComparePacking1);

    json reply;

    reply["tradeId"] = global_tradeId;

    vector<Cargo*>::iterator cargo_end = vec.end() - 1;
    if ((*cargo_end)->numofpack >= maxnbin - 100)
    {
        map<int, map<string, vector<Cargo*>>> containers;

        for (Cargo* cargo : vec){
            if(cargo->numofpack < maxnbin - 100){
                containers[cargo->numofpack][cargo->orderCode].push_back(cargo);
            }
        }

        reply["status"] = "0";
        reply["msg"] = "unenough bins";
        reply["carNum"] = containers.size();
        int vec_size = vec.size();
        for(Cargo* cargo : vec){
            if(cargo->numofpack >= maxnbin - 100)
            {
                vec_size--;
            }
        }
        reply["goodNum"] = vec_size;

        {
            int64_t totalCapacityCM3 = 0;
            double totalWeight = 0;

            for (Cargo* cargo : vec)
            {
                totalCapacityCM3 += cargo->GetVolumeCM3();
                totalWeight += cargo->GetWeight();
            }

            double totalCapacityM3 = (double)totalCapacityCM3 / 10e6;

            reply["totalCapacity"] = totalCapacityM3;
            reply["totalWeight"] = totalWeight;
        }

        reply["trainList"] = calcTrains(containers);

        return reply.dump();

    }
    else
    {
        map<int, map<string, vector<Cargo*>>> containers;

        for (Cargo* cargo : vec)
            containers[cargo->numofpack][cargo->orderCode].push_back(cargo);

        reply["status"] = "1";
        reply["carNum"] = containers.size();
        reply["goodNum"] = vec.size();

        {
            int64_t totalCapacityCM3 = 0;
            double totalWeight = 0;

            for (Cargo* cargo : vec)
            {
                totalCapacityCM3 += cargo->GetVolumeCM3();
                totalWeight += cargo->GetWeight();
            }

            double totalCapacityM3 = (double)totalCapacityCM3 / 10e6;

            reply["totalCapacity"] = totalCapacityM3;
            reply["totalWeight"] = totalWeight;
        }

        reply["trainList"] = calcTrains(containers);

        return reply.dump();
    }
}

void ShowCargo(vector<Cargo*> vec)     //װ������λ�����
{
    string str = "{\"type\" : \"bin\", \"cuboid\" : {\"x1\" : ";
    str = str + to_string((bin[0]->GetPoint1()).GetPointX());
    str = str + ", \"x2\" : ";
    str = str + to_string((bin[0]->GetPoint2()).GetPointX());
    str = str + ", \"y1\" : ";
    str = str + to_string((bin[0]->GetPoint1()).GetPointY());
    str = str + ", \"y2\" : ";
    str = str + to_string((bin[0]->GetPoint2()).GetPointY());
    str = str + ", \"z1\" : ";
    str = str + to_string((bin[0]->GetPoint1()).GetPointZ());
    str = str + ", \"z2\" : ";
    str = str + to_string((bin[0]->GetPoint2()).GetPointZ());
    str = str + "}}\n";
    // // << str;
    const char *ch = str.c_str();
    SocketUDP(ch);

    //p2_
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        Cargo* c = *it;
        // printf("��������λ��:\n");
        ShowPoint(c->GetPoint1(), c->GetPoint2());
        // << "�����ţ�" << c->GetId() << " ������: " << c->ding << " ��װ��װ�����кţ�" << c->numofpack << endl;
        //printf("{\"type\" : \"cargo\", \"cuboid\" : {\"x1\" : %d, \"x2\" : %d, \"y1\" : %d, \"y2\" : %d, \"z1\" : %d, \"z2\" : %d }}\n", (c->GetPoint1()).GetPointX(), (c->GetPoint2()).GetPointX(), (c->GetPoint1()).GetPointY(), (c->GetPoint2()).GetPointY(), (c->GetPoint1()).GetPointZ(), (c->GetPoint2()).GetPointZ());
        string str = "{\"type\" : \"cargo\", \"cuboid\" : {\"x1\" : ";
        str = str + to_string((c->GetPoint1()).GetPointX());
        str = str + ", \"x2\" : ";
        str = str + to_string((c->GetPoint2()).GetPointX());
        str = str + ", \"y1\" : ";
        str = str + to_string((c->GetPoint1()).GetPointY());
        str = str + ", \"y2\" : ";
        str = str + to_string((c->GetPoint2()).GetPointY());
        str = str + ", \"z1\" : ";
        str = str + to_string((c->GetPoint1()).GetPointZ());
        str = str + ", \"z2\" : ";
        str = str + to_string((c->GetPoint2()).GetPointZ());
        str = str + "}}\n";
        // << str;
        const char *ch = str.c_str();
        if((*it)->numofpack == CHE && ONE == 1)
        {
            //Sleep(1000);
            SocketUDP(ch);
        }
    }
}
//���������пռ�
void StateUpdate(Cargo *c, int flag)
{
    //ÿ��װ���ӽ�֮ǰ�Ŀռ�����鲢����ά��
    vece2.clear();
    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        EMSpace ems_ = *ite;
        int x1 = max((ems_.GetPoint1()).GetPointX(), (c->GetPoint1()).GetPointX());
        int y1 = max((ems_.GetPoint1()).GetPointY(), (c->GetPoint1()).GetPointY());
        int z1 = max((ems_.GetPoint1()).GetPointZ(), (c->GetPoint1()).GetPointZ());

        int x2 = min((ems_.GetPoint2()).GetPointX(), (c->GetPoint2()).GetPointX());
        int y2 = min((ems_.GetPoint2()).GetPointY(), (c->GetPoint2()).GetPointY());
        int z2 = min((ems_.GetPoint2()).GetPointZ(), (c->GetPoint2()).GetPointZ());
        //֤�������������ཻ
        if(x1 < x2 && y1 < y2 && z1 < z2)
        {
            Point p1(x1, y1, z1), p2(x2, y2, z2);
            //�����������ems(X1-X2)_, С����������p1-p2(X3-X4)
            int X1 = (ems_.GetPoint1()).GetPointX(), Y1 = (ems_.GetPoint1()).GetPointY(), Z1 = (ems_.GetPoint1()).GetPointZ();
            int X2 = (ems_.GetPoint2()).GetPointX(), Y2 = (ems_.GetPoint2()).GetPointY(), Z2 = (ems_.GetPoint2()).GetPointZ();
            int X3 = x1, Y3 = y1, Z3 = z1;
            int X4 = x2, Y4 = y2, Z4 = z2;
            //��Ϊ6���ռ�
            //1(���)
            if(X1 != X3 && Y1 != Y2 && Z1 != Z2)
            {
                Point point1(X1, Y1, Z1), point2(X3, Y2, Z2);
                int length[3];
                length[0] = X3 - X1;    //��С�ĳ���
                length[1] = Y2 - Y1;
                length[2] = Z2 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;

                    ems.max_weight = ite->max_weight;
                    ems.cargo_code = ite->cargo_code;
                    ems.layer = ite->layer;

                    vece2.push_back(ems);
                }
            }
            //2(�Ҳ�)
            if(X4 != X2 && Y1 != Y2 && Z1 != Z2)
            {
                Point point1(X4, Y1, Z1), point2(X2, Y2, Z2);
                int length[3];
                length[0] = X2 - X4;    //��С�ĳ���
                length[1] = Y2 - Y1;
                length[2] = Z2 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;

                    ems.max_weight = ite->max_weight;
                    ems.cargo_code = ite->cargo_code;
                    ems.layer = ite->layer;

                    vece2.push_back(ems);
                }
            }
            //3��ǰ�ࣩ
            if(X1 != X2 && Y1 != Y3 && Z1 != Z2)
            {
                Point point1(X1, Y1, Z1), point2(X2, Y3, Z2);
                int length[3];
                length[0] = X2 - X1;    //��С�ĳ���
                length[1] = Y3 - Y1;
                length[2] = Z2 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;

                    ems.max_weight = ite->max_weight;
                    ems.cargo_code = ite->cargo_code;
                    ems.layer = ite->layer;

                    vece2.push_back(ems);
                }
            }
            //4����ࣩ
            if(X1 != X2 && Y4 != Y2 && Z1 != Z2)
            {
                Point point1(X1, Y4, Z1), point2(X2, Y2, Z2);
                int length[3];
                length[0] = X2 - X1;    //��С�ĳ���
                length[1] = Y2 - Y4;
                length[2] = Z2 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;

                    ems.max_weight = ite->max_weight;
                    ems.cargo_code = ite->cargo_code;
                    ems.layer = ite->layer;

                    vece2.push_back(ems);
                }
            }
            //5���²ࣩ
            if(X1 != X2 && Y1 != Y2 && Z1 != Z3)
            {
                Point point1(X1, Y1, Z1), point2(X2, Y2, Z3);
                int length[3];
                length[0] = X2 - X1;    //��С�ĳ���
                length[1] = Y2 - Y1;
                length[2] = Z3 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;

                    ems.max_weight = ite->max_weight;
                    ems.cargo_code = ite->cargo_code;
                    ems.layer = ite->layer;

                    vece2.push_back(ems);
                }
            }
            //6���ϲࣩ
            if(X1 != X2 && Y1 != Y2 && Z4 != Z2)
            {
                Point point1(X1, Y1, Z4), point2(X2, Y2, Z2);
                int length[3];
                length[0] = X2 - X1;    //��С�ĳ���
                length[1] = Y2 - Y1;
                length[2] = Z2 - Z4;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (c->node[flag - 1]).l;
                    ems.yes = (c->node[flag - 1]).y;
                    ems.weight = c->weight;
                    //ems.level = (*ite).level;
                    //ems.yes = (*ite).yes;

                    if (ite->cargo_code == c->materialCode)
                    {
                        ems.cargo_code = c->materialCode;
                        ems.max_weight = ite->max_weight;
                        ems.layer = ite->layer - 1;
                        //cerr << "if - layer = " << ems.layer << endl;
                    }
                    else
                    {
                        ems.max_weight = min(ite->max_weight - c->weight, c->weight);
                        ems.cargo_code = c->materialCode;
                        //!!!						ems.layer = c->stackLevel;
                        ems.layer = (c->node[flag - 1]).stackLevel;
                        //cerr << "else - layer = " << ems.layer << endl;

                    }


                    vece2.push_back(ems);
                }
            }
        }
        //���������岻�ཻ
        else
        {
            vece2.push_back(ems_);
        }
    }
    //vece�д洢Ϊ��������Ŀռ�
    vece.clear();
    for(ite2 = vece2.begin(); ite2 != vece2.end(); ite2++)
    {
        vece.push_back(*ite2);
    }
    //��vece�д洢������ռ�ɾ����
    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        //ShowPoint(ite->GetPoint1(), ite->GetPoint2());
        for(ite2 = vece2.begin(); ite2 != vece2.end(); ite2++)
        {
            if(Compare(ite->GetPoint1(), ite2->GetPoint1()) == 2 && Compare(ite->GetPoint2(), ite2->GetPoint2()) == 2)
            {
                continue;
            }
            //ite������ite2��
            else if((Compare(ite->GetPoint1(), ite2->GetPoint1()) == 1 || Compare(ite->GetPoint1(), ite2->GetPoint1()) == 2) && (Compare(ite->GetPoint2(), ite2->GetPoint2()) == -1 || Compare(ite->GetPoint2(), ite2->GetPoint2()) == 2))
            {
                vece.erase(ite);
                ite--;
                break;
            }
        }
    }

    //��vece������Ŀռ�ɾ����
    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        Point pv_1 = ite->GetPoint1();  //�ռ������
        Point pv_2 = ite->GetPoint2();
        Point pc_1 = c->GetPoint1();    //���ӵ�����
        Point pc_2 = c->GetPoint2();
        if(pv_2.GetPointX() <= pc_1.GetPointX() && pv_1.GetPointZ() <= pc_2.GetPointZ())
        //if(pv_1.GetPointX() <= pc_1.GetPointX() && pv_1.GetPointZ() <= pc_1.GetPointZ())
        {
            vece.erase(ite);
            ite--;
        }
    }
    /*
    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        Point pv_1 = ite->GetPoint1();  //�ռ������
        Point pv_2 = ite->GetPoint2();
        Point pc_1 = c->GetPoint1();    //���ӵ�����
        Point pc_2 = c->GetPoint2();

        if(pv_1.GetPointX() + LengthOfPacking <= pc_1.GetPointX())      //�������ӷ��ù�������
        {
            //// << pv_1.GetPointX() << "    " <<  pc_1.GetPointX() << "    ";
            int x_ = pc_2.GetPointX();
            //// << x_ - LengthOfPacking << "   ";
            ite->Change(x_ - LengthOfPacking);
            //// << (ite->GetPoint1()).GetPointX() << endl;
        }
    }
    */
}

template <typename T>
T field_of(const json& object, const char* name)
{
    try
    {
        return T(object.at(name));
    }
    catch (const exception& err)
    {
        throw runtime_error(string("getting `") + name + string("`: ") + err.what());
    }
}

void JSON(string read)
{
    json input = json::parse(read);//����json����


    // string tradeId = field_of<string>(input, "tradeId"); //tradeId��
    global_tradeId = field_of<string>(input, "tradeId"); //tradeId��

    json orderList = field_of<json>(input, "orderList");//�����б�
    json binList = field_of<json>(input, "vehicleModelList");//�����б�
    int cargoid = 0;

    //// << orderList.size() << endl;
    for(int i=0; i<orderList.size(); i++) //�������ж���
    {
        json order = orderList[i];//��ǰ����

        // long long orderCode = order.at("orderCode");//����ID,string������ת
        string orderCode = field_of<string>(order, "orderCode");//����ID,string������ת

        //int orderID = stoi(tempString);
        //// << orderID << endl;
        int unloadingSequence = field_of<int>(order, "unloadingSequence");

        json boxList = field_of<json>(order,"goodList");//���Ｏ��json����
        //// << boxList.size() << endl;


        for(int j =0; j<boxList.size(); j++) //������ǰ���������л���
        {

            json box = boxList[j];//��ǰ����

            //string qty = box.at("qty"); //�������
            //int qtynum = stoi(qty);
            string materialCode = field_of<string>(box,"materialCode");//����ID,string������ת
            int qtynum = field_of<int>(box, "qty");
            //// << qtynum << endl;
            qtynum--;

            //// << idOfSet << endl;


            //string tempLength;//�õ�����ߺ�������string������ת
            /*
            tempLength=box.at("length");
            int x = stoi(tempLength);
            // << x << endl;
            tempLength=box.at("width");
            int y = stoi(tempLength);
            // << y << endl;
            tempLength=box.at("height");
            int z = stoi(tempLength);
            // << z << endl;
            tempLength=box.at("weight");
            int weight = stoi(tempLength);
            // << weight << endl;
            */
            int x = field_of<int>(box, "length");
            x = x / 10;
            int y = field_of<int>(box, "width");
            y = y/10;
            int z = field_of<int>(box, "height");
            z = z / 10;
            double weight = field_of<double>(box, "weight");

            ca[cargoid] = new Cargo(x, y, z);
            ca[cargoid]->Setnum(0);
            ca[cargoid]->SetId(cargoid);
            ca[cargoid]->ding = 10000 - unloadingSequence;
            ca[cargoid]->weight = weight;
            ca[cargoid]->materialCode = materialCode;
            ca[cargoid]->orderCode = orderCode;

            string idOfSet = "";
            if (box.contains("setCode") && !box["setCode"].is_null())
            {
                idOfSet = field_of<string>(box,"setCode");//�׻�����
                if (idOfSet != "")
                {
                    ca[cargoid]->set_code = idOfSet;
                }
            }
            //// << ca[cargoid]->weight << endl;

            json boxLimit = field_of<json>(box, "restrictionList");//װ������json����
            for(int k =0; k<boxLimit.size(); k++) //������������
            {
                json oriLimit = boxLimit[k];//��ǰ���ƹ���

                int orientation = stoi(field_of<string>(oriLimit, "flag"));//���Ʒ���
                // cerr << "input orientation = " << orientation << endl;
                orientation = map_input_orientation(orientation);
                // cerr << "map_input_orientation(orientation) = " << orientation << endl;

                int direction = 0;
                if(orientation == 1)
                    direction = 1;
                else if(orientation == 2)
                    direction = 3;
                else if(orientation == 3)
                    direction = 2;
                else if(orientation == 4)
                    direction = 4;
                else if(orientation == 5)
                    direction = 5;
                else if(orientation == 6)
                    direction = 6;

                --direction;
                //// << orientation << endl;
                //// << "     " << direction << endl;

                bool tempIsBearSurface = field_of<bool>(oriLimit,"isBear");//����������
                int isBearSurface=0;
                int bearLevel = 1;

                if(tempIsBearSurface==true)
                {
                    isBearSurface=1;
                    bearLevel = field_of<int>(oriLimit, "bearLevel");
                }


                ca[cargoid]->SetOri(direction, 0);
                (ca[cargoid]->node[direction]).l = bearLevel;
                (ca[cargoid]->node[direction]).y = isBearSurface;
                (ca[cargoid]->node[direction]).f = 1;
                ca[cargoid]->node[direction].isStack = field_of<bool>(oriLimit, "isStack");
                if (!ca[cargoid]->node[direction].isStack)
                {
                    ca[cargoid]->node[direction].stackLevel = 100;
                }
                else
                {
                    ca[cargoid]->node[direction].stackLevel = field_of<int>(oriLimit, "stackLevel");
                }
            }
            vecc.push_back(ca[cargoid]);
            cargoid++;

            while(qtynum)
            {
                //// << qtynum << endl;
                ca[cargoid] = new Cargo(x, y, z);
                ca[cargoid]->Setnum(0);
                for(int j=0; j<6; j++)
                {
                    ca[cargoid]->SetOri(j, 0);
                    (ca[cargoid]->node[j]).l = (ca[cargoid-1]->node[j]).l;
                    (ca[cargoid]->node[j]).f = (ca[cargoid-1]->node[j]).f;
                    (ca[cargoid]->node[j]).y = (ca[cargoid-1]->node[j]).y;
                    ca[cargoid]->node[j].isStack = ca[cargoid-1]->node[j].isStack;
                    ca[cargoid]->node[j].stackLevel = ca[cargoid-1]->node[j].stackLevel;
                }
                ca[cargoid]->SetId(cargoid);
                ca[cargoid]->ding = 10000 - unloadingSequence;
                ca[cargoid]->weight = weight;
                ca[cargoid]->orderCode = orderCode;
                ca[cargoid]->set_code = idOfSet;

                ca[cargoid]->materialCode = materialCode;

                vecc.push_back(ca[cargoid]);
                cargoid++;
                qtynum--;
            }

            if(x > y)
            {
                int temp = x;
                x = y;
                y = temp;
            }
            if(x > z)
            {
                int temp = x;
                x = z;
                z = temp;
            }
            if(y > z)
            {
                int temp = y;
                y = z;
                z = temp;
            }
            minCargo1 = min(minCargo1, x);
            minCargo2 = min(minCargo2, y);
            minCargo3 = min(minCargo3, z);
        }
    }


    int binid = 0;

    for(int i=0; i<binList.size(); i++) //���������б�
    {

        json binq = binList[i];//��ǰ����

        //json binAttribute = bin.at("��������");//��������

        //string idOfBin = binAttribute.at("�����ͺ�");//�����ͺ�
        string typeOfBin = field_of<string>(binq,"modelCode");//��������

        /*
        string intTemp=bin.at("weight");//תint��ʱ����
        int weightOfBin = stoi(intTemp);//��������

        intTemp = bin.at("length");
        int x = stoi(intTemp);//���䳤

        intTemp = bin.at("width");
        int y = stoi(intTemp);//�����

        intTemp = bin.at("height");
        int z = stoi(intTemp);//�����

        string binNumbers = bin.at("qty");//�����������
        int numOfBin = stoi(binNumbers);
        */
        int weightOfBin = field_of<int>(binq, "weight");
        int x = field_of<int>(binq, "length");
        x = x / 10;
        int y = field_of<int>(binq, "width");
        y = y / 10;
        int z = field_of<int>(binq,"height");
        z = z / 10;
        int numOfBin = field_of<int>(binq, "qty");

        while(numOfBin)
        {
            Point p1(0, 0, 0);
            Point p2(x, y, z);
            bin[binid] = new Bin(p1, p2);
            //bin[binid] = new Bin(p1, p2);
            bin[binid]->SetNum(binid);
            bin[binid]->modelCode = typeOfBin;
            binid++;
            numOfBin--;
        }
    }
    int nbin = binid;
    for (int i = 0; i < nbin; i++)
    {
        for (int j = i + 1; j < nbin; j++)
        {
            if (bin[i]->GetHigh() < bin[j]->GetHigh())
            {
                Point p1 = bin[i]->GetPoint1();
                Point p2 = bin[i]->GetPoint2();
                string strpp = bin[i]->modelCode;

                bin[i]->SetPoint(bin[j]->GetPoint1(), bin[j]->GetPoint2());
                bin[i]->modelCode = bin[j]->modelCode;

                bin[j]->SetPoint(p1, p2);
                bin[j]->modelCode = strpp;

            }
        }
    }
    for (int i = 0; i < nbin; i++)
    {
        bin[i]->SetNum(i);
    }

    //�¼�100���ֿ�
    for(int i = 0; i<100; i++)
    {
        Point max_p1(0, 0, 0);
        Point max_p2(1000, 1000, 1000);
        bin[nbin] = new Bin(max_p1, max_p2);
        bin[nbin]->SetNum(nbin);
        bin[nbin]->modelCode = "storehouse";
        nbin++;
    }


    maxn = cargoid;
    maxnbin = nbin;
}

string binPackingAlgorithm(const char* str)
{
    Init();
    JSON(str);

    Chromosome* chromosome_ans = GeneticAlgorithm();
    bool valid = chromosome_ans && chromosome_ans->PackingOf(vecc);

    string ans = another_ANS(vecc, valid);

    // char* reply = (char*)malloc(ans.size());
    // memcpy(reply, ans.data(), ans.size());
    // strcpy(reply, ans.c_str());

    delete_vector(bin, maxnbin);

    delete_vector(ca, maxn);

    for (auto ptr : allocated_chrosomes)
        delete ptr;

    // return reply;
    return ans;
}

int main()
{
    Init();

    ostringstream os;
    os << cin.rdbuf();

    cout << binPackingAlgorithm(os.str().c_str()) << endl;
}
