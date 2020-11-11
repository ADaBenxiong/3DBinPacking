#ifndef CHROMOSOME_H_INCLUDED
#define CHROMOSOME_H_INCLUDED

#include"Cargo.h"
// #include <Winsock2.h>
// #include <windows.h>
// #include<bits/stdc++.h>
#include <vector>
//#include"stdafx.h"
using namespace std;
static int maxn = 80;    //表示货物的数量
static const int maxn_ = 5000;
static int maxnbin = 10;  //表示集装箱的数量
static const int maxnbin_ = 1000;
static int ID = 0;

class Chromosome
{
public:
    Chromosome(double *variable)
		: ch(maxn_ * 2)
    {
        for(int i=0; i<maxn * 2; i++)
        {
            ch[i] = variable[i];
        }
        this->fitness = -1;
        this->refitness = 0;
        this->sumfitness = 0;
    }
    ~Chromosome(){}
    void SetValue(int i, double k);    //染色体基因赋值
    double GetValue(int i) const;       //得到染色体基因值

    double Fitness(vector<Cargo*> vec, Bin *bin);   //计算适应值
    bool PackingOf(vector<Cargo*> vec);   //根据染色体值装箱显示

    void SetId(int k);   //染色体编号赋值
    int GetId() const;  //获取染色体编号
    void SetFitness(double k);      //适应值赋值操作
    double GetFitness() const;
    void SetReFitness(double k);    //计算适应概率
    double GetReFitness() const;    //得到染色体适应概率
    void SetSumFitness(double k);   //（轮盘赌算法）计算适应概率总和
    double GetSumFitness() const;  //得到染色体适应概率总和
    int demo;   //表示这件物品订单号
    int level;  //表示物品的称重级别
private:
    int id; //染色体的编号
    // double ch[maxn_ * 2]; //表示染色体的基因值
	vector<double> ch;
    double fitness = -1; //表示该条染色体的适应值
    double refitness;   //表示染色体的适应概率
    double sumfitness;  //表示累加的概率
};

void Chromosome::SetId(int k)
{
    id = k;
}
int Chromosome::GetId() const
{
    return id;
}
void Chromosome::SetValue(int i, double k)
{
    ch[i] = k;
}
double Chromosome::GetValue(int i) const
{
    return ch[i];
}
void Chromosome::SetFitness(double k)
{
    fitness = k;
}
double Chromosome::GetFitness() const
{
    return fitness;
}
void Chromosome::SetReFitness(double k)
{
    refitness = k;
}
double Chromosome::GetReFitness() const
{
    return refitness;
}
void Chromosome::SetSumFitness(double k)
{
    sumfitness = k;
}
double Chromosome::GetSumFitness() const
{
    return sumfitness;
}

#endif // CHROMOSOME_H_INCLUDED
