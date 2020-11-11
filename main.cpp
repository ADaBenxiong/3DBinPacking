

//问题1：虽然对套机进行了捆绑操作，但是不确定是否会出现恰好装于两车的情况
//问题2：悬空问题通过不断的迭代来达到更好的效果，但是也有可能出现悬空的问题
//问题3：多箱问题对于如何选择集装箱并未做出考虑，需要通过经济效益判断如何选箱
//问题4：给出的真实数据存在误差，有问题（高度过高的冰箱，以及装箱率问题）

#include <jni.h>
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

template <typename Container>
void delete_vector(const Container& container, int size) {
	for (int i = 0; i < size; ++i)
		delete container[i];
}

template <typename Container>
void delete_vector(const Container& container) {
	for (auto ptr : container)
		delete ptr;
}

using namespace std;
using json = nlohmann::json;
//链接静态库
#pragma comment( lib, "ws2_32.lib" )
/*
UDP：
初始化网络库
1.创建套接字
2.绑定套接字
3.接收数据
4.关闭网络
*/
//DWORD WINAPI RecvFun(LPVOID lp);
int SUM = 0;        //表示遗传算法迭代次数
int CHE = 0;    //表示显示车的编号
int ONE = 1;    //表示一个一个显示
int minCargo1 = 2000; //剪掉空间过小的ems(分别表示第一第二第三小)
int minCargo2 = 2000;
int minCargo3 = 2000;

//================================================集装箱以及货物初始信息
Point p1_(0, 0, 0);
Point p2_(720, 230, 215);

string global_tradeId;

Bin* bin[maxnbin_]; //集装箱的数量
Cargo *ca[maxn_];   //货物的数量
//double SUM[maxnbin_];   //每个车的装箱率
int LengthOfPacking = 100;
/*
void init()
{
//集装箱的产生
    for(int i = 0; i < maxnbin; i++)
    {
        bin[i] = new Bin(p1_, p2_);
        bin[i]->SetNum(i);
    }
}
*/



//货物的信息存储
vector<Cargo*>vecc;
vector<Cargo*>::iterator itc;
//最大空闲空间序列
vector<EMSpace>vece;
//用于维护最大空闲空间
vector<EMSpace>vece2;
vector<EMSpace>::iterator ite;
vector<EMSpace>::iterator ite2;

//================================================遗传算法函数以及变量
const int MAXN = 20;    //表示遗传算法每一代种群的数量
const int EV = 0; //表示遗传算法进化的代数
const double SaProbability = 0.100; //优秀染色体保留的概率
const double OvProbability = 0.600; //交叉的概率
const double VaProbability = 0.300;  //变异的概率
const double OvElite = 0.500;   //杂交过程选择精英概率
const double VaGene = 0.400;    //发生染色体变异的概率

vector<Chromosome*>vec_now;    //进化种群
vector<Chromosome*>vec_next;    //进化下一代种群
vector<Chromosome*>vec_elite;   //精英染色体种群
Chromosome* chrom[MAXN];        //染色体数组（最初的生成）

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

void Initialize();  //随机初始化种群， 得到第一代个体
void CaculaFitness();   //计算个体的适应值

void Selecting();     //种群选择
void Crossing();    //杂交
void Variating();   //变异
Chromosome* GeneticAlgorithm();    //遗传算法(最后返回一条最佳的染色体)

//==============================================装箱算法函数以及变量

void ShowPoint(Point p1, Point p2);     //显示箱体的两个顶点
int Compare(Point p1, Point p2);        //比较两个点的位置， 返回值1表示p1在p2的右后上方
bool  DFTRC_2(vector<Cargo*>vec, Cargo *c, Bin *b, EMSpace*);      //箱子的编号、装箱空间的选择
// EMSpace DFTRC_2(vector<Cargo*>vec, Cargo *c, Bin *b);      //箱子的编号、装箱空间的选择

int CargoOrientation(vector<Cargo*>vec, Cargo *c, EMSpace ems);    //装箱方向的选择
void CargoPacking(Cargo *c, EMSpace ems, int flag);   //用于装箱的箱子， 装箱的空间， 箱子的方向向量元素
void StateUpdate(Cargo *c, int flag);         //维护装箱的空间
void ShowCargo(vector<Cargo*> vec);     //装箱箱子位置输出
//void SocketUDP(const char* sendBuf);       //通过UDP进行通信
//发送UDP通信过程
void SocketUDP(const char* sendBuf)
{
    // WSADATA wsaData;
    // //初始化socket
    // WSAStartup(MAKEWORD(2,2),&wsaData);
    // //套接字初始化
    // SOCKET sendSocket;
    // sendSocket=socket(AF_INET,SOCK_DGRAM,IPPROTO_UDP);
    // //设置即将连服务器地址
    // sockaddr_in seAddr;
    // seAddr.sin_family=AF_INET;
    // seAddr.sin_port=htons(55554);
    // seAddr.sin_addr.S_un.S_addr=inet_addr("127.0.0.1");//;htonl(INADDR_ANY)

    // //初始化
    // int bufLen=strlen(sendBuf);
    // //绑定
    // sendto(sendSocket,sendBuf,bufLen,0,(SOCKADDR *)&seAddr,sizeof(seAddr));
    // //发送完成，关闭socket
    // closesocket(sendSocket);
    // //释放资源并退出
    // WSACleanup();
}
//===============================================遗传算法类的方法的实现
bool ComparePacking1(Cargo* x, Cargo* y)     //将箱子按照染色体基因进行排序
{
    if(x->GetBPS() != y->GetBPS())
        return x->GetBPS() < y->GetBPS();
    else
        return x->GetId() < y->GetId();
}
bool ComparePacking2(Cargo* x, Cargo* y)    //将箱子按照序号进行排序
{
    return x->GetId() < y->GetId();
}
double Chromosome::Fitness(vector<Cargo*> vec, Bin* b)        //染色体的适应值计算函数
{
    //集装箱初始信息维护
    vece.clear();
    EMSpace ems(b->GetPoint1(), b->GetPoint2());
    vece.push_back(ems);
    sort(vec.begin(), vec.end(), ComparePacking2);

    int k = 0;  //表示染色体基因值序列
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        (*it)->SetBPS(ch[k] + (*it)->ding);
        (*it)->SetVBO(ch[maxn + k]);
        k++;
    }

    sort(vec.begin(), vec.end(), ComparePacking1);

    int numk = 0;   //表示集装箱的编号
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        //cout << (*it)->GetLength() << " " << (*it)->GetWidth() << " " << (*it)->GetHigh() << endl;
        // EMSpace ems_select = DFTRC_2(vec, *it, bin[numk]);   //选择用于装箱的空间
		EMSpace ems_select(ems);
        bool valid = DFTRC_2(vec, *it, bin[numk], &ems_select);   //选择用于装箱的空间
		if (!valid)
			return 0;

        //cout << ems_select.GetLength() << " " << ems_select.GetWidth() << " " << ems_select.GetHigh() << endl;
        numk = (*it)->numofpack;
        //cout << numk;
        int flag = CargoOrientation(vec, *it, ems_select); //返回即为箱子的方向
        //cout << flag << endl;
        CargoPacking(*it, ems_select, flag);  //进行装箱
        StateUpdate(*it, flag);
    }
    vece.clear();//装箱空间最后清除
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
    double value_k = numk + (((max_x) * 1.0) / ((bin[numk]->GetPoint2()).GetPointX())* 1.0);  //是一个浮点数，表示装货一共需要多少个箱子
    //double value_k = numk + (((max_x * max_y * max_z) * 1.0) / (p2_.GetPointX() * p2_.GetPointY() * p2_.GetPointZ())* 1.0);  //是一个浮点数，表示装货一共需要多少个箱子
    value_k = 1.0 / value_k;
    fitness = value_k;
    sort(vec.begin(), vec.end(), ComparePacking2);
    return value_k;
}

bool Chromosome::PackingOf(vector<Cargo*> vec)
{

    printf("装箱显示\n");
    //集装箱初始信息维护
    vece.clear();
    EMSpace ems(bin[0]->GetPoint1(), bin[0]->GetPoint2());
    vece.push_back(ems);
    sort(vec.begin(), vec.end(), ComparePacking2);
    int k = 0;  //表示染色体基因值序列
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
		int ding = (*it)->ding;
        (*it)->SetBPS(ch[k] + ding);
        (*it)->SetVBO(ch[maxn + k]);
        k++;
    }
    sort(vec.begin(), vec.end(), ComparePacking1);

    int numk = 0;   //表示集装箱的编号                                          //表示装载于哪个集装箱
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        // EMSpace ems_select = DFTRC_2(vec, *it, bin[numk]);   //选择用于装箱的空间
		EMSpace ems_select(ems);
        bool valid = DFTRC_2(vec, *it, bin[numk], &ems_select);   //选择用于装箱的空间
		if (!valid)
			return false;

        numk = (*it)->numofpack;
        cout << numk;
        int flag = CargoOrientation(vec, *it, ems_select); //返回即为箱子的方向

        (*it)->flag = flag;
        if((*it)->flag == 2)
            (*it)->flag = 3;
        else if((*it)->flag == 3)
            (*it)->flag = 2;
        //cout << flag << endl;
        CargoPacking(*it, ems_select, flag);  //进行装箱
        StateUpdate(*it, flag);
        //sum_all += (*it)->GetLength() * (*it)->GetWidth() * (*it)->GetHigh();
    }
    vece.clear();//装箱空间最后清除
    //将装箱的顺序和位置返回

    sort(vec.begin(), vec.end(), ComparePacking2);

    double SUM[maxnbin_];
    memset(SUM, 0, sizeof(SUM));
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        int np = (*it)->numofpack;
        SUM[np] += (*it)->GetLength() * (*it)->GetWidth() * (*it)->GetHigh();
    }
    for(int i=0; i<=numk; i++)
    {
        cout << "集装箱" << i << "装箱率：" << SUM[i] * 1.0 / (bin[i]->GetLength() * bin[i]->GetWidth() * bin[i]->GetHigh()) * 1.0 << endl;
    }
    // ShowCargo(vec);

	return true;
}
//=======================================================遗传算法函数的实现

void Initialize()   //随机初始化种群， 得到第一代种群
{
    vec_now.clear();
    srand(time(NULL));  //产生随机的种子
    ID = 0;
    for(int i=0; i<MAXN; i++)
    {
        double ch_variable[maxn * 2];
        for(int j=0; j<maxn * 2; j++)
        {
            double p = (rand() % 1000)/ 1000.0;
            ch_variable[j] = p;
        }/*
        for(int j = 0; j<maxn * 2; j++)
        {
            cout << ch_variable[j] << " ";
        }
        cout << endl;
        */
        chrom[i] = new Chromosome(ch_variable);
		allocated_chrosomes.push_back(chrom[i]);
        chrom[i]->SetId(ID);
        ID++;
        vec_now.push_back(chrom[i]);
    }
    cout << "染色体容器大小:" << vec_now.size() << endl;
}
void CaculaFitness()    //计算适应值概率每一代所有染色体都被计算
{
    double sum = 0; //适应值概率总和
    double value[vec_now.size()];
    int k = 0;  //用于对value值进行遍历
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)
    {
        //if((*it)->GetFitness() == -1)
        value[k] = (*it)->Fitness(vecc, bin[0]);    //从第0个箱子开始
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
        (*it)->SetReFitness(m);     //概率越大越好
        //cout << 1.0 / (*it)->GetFitness() * 1.0 << endl;
        max_ = max((*it)->GetFitness(), max_);                                                                 //用于调试显示需要迭代足够多的次数
        k++;
    }
    SUM++;
    cout << SUM << " " << max_ << endl;
}
bool CompareChromosome(Chromosome* x, Chromosome* y)    //染色体排序算法
{
    if(x->GetReFitness() != y->GetReFitness())
        return x->GetReFitness() > y->GetReFitness();
    else
        return x->GetId() < y->GetId();
}
//选择精英染色体保留
void Selecting()
{
    vec_next.clear();
    vec_elite.clear();

    sort(vec_now.begin(), vec_now.end(), CompareChromosome);
    int Length = (int)(vec_now.size() * SaProbability);    //精英染色体的个数
    int k = 0;
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it != vec_now.end(); it++)
    {
        if(k < Length)
        {
            vec_next.push_back(*it);        //将精英染色体选择到下一代中
            vec_elite.push_back(*it);
            //cout << "精英" << (*it)->GetReFitness() << endl << endl;
        }
        else
            break;
        k++;
    }
    //轮盘赌算法 赋予染色体概率值
    double sum = 0;
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)
    {
        sum = sum + (*it)->GetReFitness();
        (*it)->SetSumFitness(sum);
        //cout <<(*it)->GetSumFitness() << endl;
    }
    //cout << vec_next.size() << endl;
}
//杂交产生染色体
void Crossing()
{
    //cout << "Test2" << endl;
    srand(time(NULL));
    int Length = vec_elite.size();    //精英染色体的个数
    for(int j = 0; j < (int)(vec_now.size() * OvProbability); j++)  //杂交的过程
    {
        Chromosome* chrom1 = vec_now.front(), *chrom2 = vec_now.front();
        int k_ans = rand()%Length; //被选择的精英染色体的编号
        int k_now = 0;  //表示正在遍历染色体的编号
        for(vector<Chromosome*>::iterator it = vec_elite.begin(); it!= vec_elite.end(); it++)//选择第一条染色体
        {
            if(k_now == k_ans)
            {
                chrom1 = *it;
                break;
            }
            k_now++;
        }
        double value = (rand()%1000 + 1)/1000.0;    //轮盘赌选择的概率值
        for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)       //选择第二条染色体
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
            //cout << chrom1->GetReFitness() << " " << chrom2->GetReFitness() << endl;
            double ans[maxn * 2];
            for(int i=0; i<maxn * 2; i++)
            {
                value = (rand()%1000 + 1)/1000.0;    //轮盘赌选择的概率值
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
    //cout << vec_next.size() << endl;
}
//变异产生染色体
void Variating()
{
    //cout << "Test3" << endl;
    srand(time(NULL));
    int Length = vec_now.size() - vec_next.size();  //变异染色体的数量
    //cout << Length << endl;
    for(int i=0; i<Length; i++)
    {
        Chromosome* chrom = vec_now.front();   //选择变异的染色体
        double value = (rand()%1000 + 1)/1000.0;    //轮盘赌选择的概率值
        for(vector<Chromosome*>::iterator it = vec_now.begin(); it!= vec_now.end(); it++)       //选择变异染色体
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
            value = (rand()%1000 + 1)/1000.0;    //轮盘赌选择的概率值
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
    //cout << vec_now.size() << endl;
}
Chromosome* GeneticAlgorithm() //遗传算法
{
    Initialize();   //初始化种群，随机生成第一代个体

    for(int i = 0; i < EV; i++) //每一代进行遗传算法
    {
        CaculaFitness();    //计算每条染色体适应值概率
        Selecting();    //选择染色体
        Crossing();     //杂交染色体
        Variating();    //变异染色体
    }

    CaculaFitness();   //计算最后一代染色体的适应值概率

    double max_ = 0;
    Chromosome* ch = nullptr; //最终被选择的染色体
    for(vector<Chromosome*>::iterator it = vec_now.begin(); it != vec_now.end(); it++)
    {
        //cout << (*it)->GetFitness() << endl;
        if((*it)->GetReFitness() > max_)
        {
            max_ = (*it)->GetReFitness();
            ch = *it;
        }
    }
    return ch;
}
//==========================================================装箱程序
//将箱子两点显示
void ShowPoint(Point p1, Point p2)
{
    cout << p1.GetPointX() << " " << p1.GetPointY() << " " << p1.GetPointZ() << endl;
    cout << p2.GetPointX() << " " << p2.GetPointY() << " " << p2.GetPointZ() << endl;
}
int Compare(Point p1, Point p2)
{
    int k = 0;
    if(p1.GetPointX() >= p2.GetPointX() && p1.GetPointY() >= p2.GetPointY() && p1.GetPointZ() >= p2.GetPointZ())
    {
        k = 1;   //说明点1在点2的右后上侧
    }
    if (p1.GetPointX() <= p2.GetPointX() && p1.GetPointY() <= p2.GetPointY() && p1.GetPointZ() <= p2.GetPointZ())
    {
        k = -1;   //说明点2在点1的右后上侧
    }
    if(p1.GetPointX() == p2.GetPointX() && p1.GetPointY() == p2.GetPointY() && p1.GetPointZ() == p2.GetPointZ())
    {
        k = 2;   //说明点1和点2相同
    }
    return k;
}

//选择用于装箱的空间(返回数据仍存在问题)
// EMSpace DFTRC_2(vector<Cargo*>vec, Cargo *c, Bin *b)
bool DFTRC_2(vector<Cargo*>vec, Cargo *c, Bin *b, EMSpace* space)
{
    //cout << endl;
    int k = b->GetNum();    // 箱子编号
    double maxdist = -1;
    vector<EMSpace>::iterator it_ems = vece.end();
    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        //cout <<endl;
        int x = (ite->GetPoint1()).GetPointX();
        int y = (ite->GetPoint1()).GetPointY();
        int z = (ite->GetPoint1()).GetPointZ();
        //方向1
        //if(ite->GetLength() >= c->GetLength() && ite->GetWidth() >= c->GetWidth() && ite->GetHigh() >= c->GetHigh())
        if((c->node[0]).f == 1 &&c->weight <= ite->weight && ite->GetLength() >= c->GetLength() && ite->GetWidth() >= c->GetWidth() && ite->GetHigh() >= c->GetHigh() && (c->node[0]).l >= ite->level && ite->yes == 1)
        {
            int k = 0;  //若k = 5表示可以放得下
            Point point1(x, y, z);
            Point point2(x + c->GetLength(), y + c->GetWidth(), z + c->GetHigh());  //装箱的两点
            if(z == 0)
            {
                k = 5;
            }
            else    //堆积问题
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
                //cout << "        1" << endl;
            }
        }

        //方向2
        //if(ite->GetLength() >= c->GetLength() && ite->GetWidth() >= c->GetHigh() && ite->GetHigh() >= c->GetWidth())
        if((c->node[1]).f == 1 &&c->weight <= ite->weight&& ite->GetLength() >= c->GetLength() && ite->GetWidth() >= c->GetHigh() && ite->GetHigh() >= c->GetWidth() && (c->node[1]).l >= ite->level && ite->yes == 1)
        {
            int k = 0;  //如k = 5表示可以放得下
            Point point1(x, y, z);
            Point point2(x + c->GetLength(), y + c->GetHigh(), z + c->GetWidth());  //装箱的两点
            if(z == 0)
            {
                k = 5;
            }
            else    //堆积问题
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
                //cout << point[0]->GetPointX() << " " << point[0]->GetPointY() << " " << point[0]->GetPointZ() << endl;
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
                //cout << "        2" << endl;
            }
        }

        //方向3
        //if(ite->GetLength() >= c->GetWidth() && ite->GetWidth() >= c->GetLength() && ite->GetHigh() >= c->GetHigh())
        if((c->node[2]).f == 1  &&c->weight <= ite->weight && ite->GetLength() >= c->GetWidth() && ite->GetWidth() >= c->GetLength() && ite->GetHigh() >= c->GetHigh() && (c->node[2]).l >= ite->level && ite->yes == 1)
        {
            int k = 0;  //如k = 5表示可以放得下
            Point point1(x, y, z);
            Point point2(x + c->GetWidth(), y + c->GetLength(), z + c->GetHigh());  //装箱的两点
            if(z == 0)
            {
                k = 5;
            }
            else    //堆积问题
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
                //cout << "        3" << endl;
            }
        }
        //方向4
        //if(ite->GetLength() >= c->GetWidth() && ite->GetWidth() >= c->GetHigh() && ite->GetHigh() >= c->GetLength())
        if((c->node[3]).f == 1 &&c->weight <= ite->weight && ite->GetLength() >= c->GetWidth() && ite->GetWidth() >= c->GetHigh() && ite->GetHigh() >= c->GetLength() && (c->node[3]).l >= ite->level && ite->yes == 1)
        {
            int k = 0;  //如k = 5表示可以放得下
            Point point1(x, y, z);
            Point point2(x + c->GetWidth(), y + c->GetHigh(), z + c->GetLength());  //装箱的两点
            if(z == 0)
            {
                k = 5;
            }
            else    //堆积问题
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
                //cout << "        4" << endl;
            }
        }
        //方向5
        //if(ite->GetLength() >= c->GetHigh() && ite->GetWidth() >= c->GetLength() && ite->GetHigh() >= c->GetWidth())
        if((c->node[4]).f == 1 &&c->weight <= ite->weight && ite->GetLength() >= c->GetHigh() && ite->GetWidth() >= c->GetLength() && ite->GetHigh() >= c->GetWidth() && (c->node[4]).l >= ite->level && ite->yes == 1)
        {
            int k = 0;  //如k = 5表示可以放得下
            Point point1(x, y, z);
            Point point2(x + c->GetHigh(), y + c->GetLength(), z + c->GetWidth());  //装箱的两点
            if(z == 0)
            {
                k = 5;
            }
            else    //堆积问题
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
        //方向6
        //if(ite->GetLength() >= c->GetHigh() && ite->GetWidth() >= c->GetWidth() && ite->GetHigh() >= c->GetLength())
        if((c->node[5]).f == 1 &&c->weight <= ite->weight && ite->GetLength() >= c->GetHigh() && ite->GetWidth() >= c->GetWidth() && ite->GetHigh() >= c->GetLength() && (c->node[5]).l >= ite->level && ite->yes == 1)
        {
            int k = 0;  //如k = 5表示可以放得下
            Point point1(x, y, z);
            Point point2(x + c->GetHigh(), y + c->GetWidth(), z + c->GetLength());  //装箱的两点
            if(z == 0)
            {
                k = 5;
            }
            else    //堆积问题
            {
                Point *point[5];
                point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
                point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
                point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
                int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
        // return *it_ems; //返回用于装箱的空间
		*space = *it_ems;
		return true;
    }
}

//确定装箱的方向
int CargoOrientation(vector<Cargo*>vec, Cargo *c, EMSpace ems)
{
    c->Setnum(0);
    for(int j=0; j<6; j++)
    {
        c->SetOri(j, 0);
    }
    //cout << (c->node[0]).f << " " << c->weight << " " << ems.weight << " " << (c->node[0]).l << " " << ems.level << " " << ems.yes << endl;
    //方向1
    //if(ems.GetLength() >= c->GetLength() && ems.GetWidth() >= c->GetWidth() && ems.GetHigh() >= c->GetHigh())
    if((c->node[0]).f == 1&& c->weight <= ems.weight && (c->node[0]).l >= ems.level && ems.yes == 1 && ems.GetLength() >= c->GetLength() && ems.GetWidth() >= c->GetWidth() && ems.GetHigh() >= c->GetHigh())
    {
        //cout << "yes" << endl;
        int k = 0;  //若k = 5表示可以放得下
        int x = (ems.GetPoint1()).GetPointX();
        int y = (ems.GetPoint1()).GetPointY();
        int z = (ems.GetPoint1()).GetPointZ();

        Point point1(x, y, z);
        Point point2(x + c->GetLength(), y + c->GetWidth(), z + c->GetHigh());  //装箱的两点
        if(z == 0)
        {
            k = 5;
        }
        else    //堆积问题
        {
            Point *point[5];
            point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
            point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
            int k_ = c->Getnum();   //可放箱子个数选择
            c->SetOri(k_, 1);
            c->Setnum(k_ + 1);
        }
    }

    //方向2
    //if(ems.GetLength() >= c->GetLength() && ems.GetWidth() >= c->GetHigh() && ems.GetHigh() >= c->GetWidth())
    if((c->node[1]).f == 1 && c->weight <= ems.weight &&(c->node[1]).l >= ems.level && ems.yes == 1&& ems.GetLength() >= c->GetLength() && ems.GetWidth() >= c->GetHigh() && ems.GetHigh() >= c->GetWidth())
    {
        int k = 0;  //若k = 5表示可以放得下
        int x = (ems.GetPoint1()).GetPointX();
        int y = (ems.GetPoint1()).GetPointY();
        int z = (ems.GetPoint1()).GetPointZ();

        Point point1(x, y, z);
        Point point2(x + c->GetLength(), y + c->GetHigh(), z + c->GetWidth());  //装箱的两点
        if(z == 0)
        {
            k = 5;
        }
        else    //堆积问题
        {
            Point *point[5];
            point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
            point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
            int k_ = c->Getnum();   //可放箱子个数选择
            c->SetOri(k_, 2);
            c->Setnum(k_ + 1);
        }
    }

    //方向3
    //if(ems.GetLength() >= c->GetWidth() && ems.GetWidth() >= c->GetLength() && ems.GetHigh() >= c->GetHigh())
    if((c->node[2]).f == 1 && c->weight <= ems.weight &&(c->node[2]).l >= ems.level && ems.yes == 1&& ems.GetLength() >= c->GetWidth() && ems.GetWidth() >= c->GetLength() && ems.GetHigh() >= c->GetHigh())
    {
        int k = 0;  //若k = 5表示可以放得下
        int x = (ems.GetPoint1()).GetPointX();
        int y = (ems.GetPoint1()).GetPointY();
        int z = (ems.GetPoint1()).GetPointZ();

        Point point1(x, y, z);
        Point point2(x + c->GetWidth(), y + c->GetLength(), z + c->GetHigh());  //装箱的两点
        if(z == 0)
        {
            k = 5;
        }
        else    //堆积问题
        {
            Point *point[5];
            point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
            point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
            int k_ = c->Getnum();   //可放箱子个数选择
            c->SetOri(k_, 3);
            c->Setnum(k_ + 1);
        }
    }
    //方向4
    //if(ems.GetLength() >= c->GetWidth() && ems.GetWidth() >= c->GetHigh() && ems.GetHigh() >= c->GetLength())
    if((c->node[3]).f == 1 && c->weight <= ems.weight &&(c->node[3]).l >= ems.level && ems.yes == 1 && ems.GetLength() >= c->GetWidth() && ems.GetWidth() >= c->GetHigh() && ems.GetHigh() >= c->GetLength())
    {
        int k = 0;  //若k = 5表示可以放得下
        int x = (ems.GetPoint1()).GetPointX();
        int y = (ems.GetPoint1()).GetPointY();
        int z = (ems.GetPoint1()).GetPointZ();

        Point point1(x, y, z);
        Point point2(x + c->GetWidth(), y + c->GetHigh(), z + c->GetLength());  //装箱的两点
        if(z == 0)
        {
            k = 5;
        }
        else    //堆积问题
        {
            Point *point[5];
            point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
            point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
            int k_ = c->Getnum();   //可放箱子个数选择
            c->SetOri(k_, 4);
            c->Setnum(k_ + 1);
        }
    }
    //方向5
    //if(ems.GetLength() >= c->GetHigh() && ems.GetWidth() >= c->GetLength() && ems.GetHigh() >= c->GetWidth())
    if((c->node[4]).f == 1 && c->weight <= ems.weight &&(c->node[4]).l >= ems.level && ems.yes == 1 && ems.GetLength() >= c->GetHigh() && ems.GetWidth() >= c->GetLength() && ems.GetHigh() >= c->GetWidth())
    {
        int k = 0;  //若k = 5表示可以放得下
        int x = (ems.GetPoint1()).GetPointX();
        int y = (ems.GetPoint1()).GetPointY();
        int z = (ems.GetPoint1()).GetPointZ();

        Point point1(x, y, z);
        Point point2(x + c->GetHigh(), y + c->GetLength(), z + c->GetWidth());  //装箱的两点
        if(z == 0)
        {
            k = 5;
        }
        else    //堆积问题
        {
            Point *point[5];
            point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
            point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
            int k_ = c->Getnum();   //可放箱子个数选择
            c->SetOri(k_, 5);
            c->Setnum(k_ + 1);
        }
    }
    //方向6
    //if(ems.GetLength() >= c->GetHigh() && ems.GetWidth() >= c->GetWidth() && ems.GetHigh() >= c->GetLength())
    if((c->node[5]).f == 1 && c->weight <= ems.weight &&(c->node[5]).l >= ems.level && ems.yes == 1 && ems.GetLength() >= c->GetHigh() && ems.GetWidth() >= c->GetWidth() && ems.GetHigh() >= c->GetLength())
    {
        int k = 0;  //若k = 5表示可以放得下
        int x = (ems.GetPoint1()).GetPointX();
        int y = (ems.GetPoint1()).GetPointY();
        int z = (ems.GetPoint1()).GetPointZ();

        Point point1(x, y, z);
        Point point2(x + c->GetHigh(), y + c->GetWidth(), z + c->GetLength());  //装箱的两点
        if(z == 0)
        {
            k = 5;
        }
        else    //堆积问题
        {
            Point *point[5];
            point[0] = new Point((point1.GetPointX() + point2.GetPointX())/2, (point1.GetPointY() + point2.GetPointY())/2, point1.GetPointZ()); //五点判断法
            point[1] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[2] = new Point((3*point1.GetPointX() + point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            point[3] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (3*point1.GetPointY() + point2.GetPointY())/4, point1.GetPointZ());
            point[4] = new Point((point1.GetPointX() + 3*point2.GetPointX())/4, (point1.GetPointY() + 3*point2.GetPointY())/4, point1.GetPointZ());
            int a[5] = {0, 0, 0, 0, 0}; //表示五点均未支撑

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
            int k_ = c->Getnum();   //可放箱子个数选择
            c->SetOri(k_, 6);
            c->Setnum(k_ + 1);
        }
    }

    int flag = (int)(c->GetVBO() * c->Getnum());
    return c->GetOri(flag);
}

//进行装箱操作
void CargoPacking(Cargo *c, EMSpace ems, int flag)   //用于装箱的箱子， 装箱的空间， 箱子的方向向量元素
{
    Point cargop1 = ems.GetPoint1();
    c->SetPoint1(cargop1);
    int x_ = (ems.GetPoint1()).GetPointX();
    int y_ = (ems.GetPoint1()).GetPointY();
    int z_ = (ems.GetPoint1()).GetPointZ();

    //cout << "test 1" << endl;
    //cout << x_ << " " << y_ << " " << z_ << endl;
    //cout << "test 2" << endl;
    //cout << c->GetLength() << " " << c->GetWidth() << " " << c->GetHigh() << endl;
    //cout << "test 3" << endl;
    if(flag == 1)
    {
        int x = x_ + c->GetLength();
        int y = y_ + c->GetWidth();
        int z = z_ + c->GetHigh();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //cout << x << " " << y << " " << z << endl;
    }
    else if(flag == 2)
    {
        int x = x_ + c->GetLength();
        int y = y_ + c->GetHigh();
        int z = z_ + c->GetWidth();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //cout << x << " " << y << " " << z << endl;
    }
    else if(flag == 3)
    {
        int x = x_ + c->GetWidth();
        int y = y_ + c->GetLength();
        int z = z_ + c->GetHigh();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //cout << x << " " << y << " " << z << endl;
    }
    else if(flag == 4)
    {
        int x = x_ + c->GetWidth();
        int y = y_ + c->GetHigh();
        int z = z_ + c->GetLength();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //cout << x << " " << y << " " << z << endl;
    }
    else if(flag == 5)
    {
        int x = x_ + c->GetHigh();
        int y = y_ + c->GetLength();
        int z = z_ + c->GetWidth();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //cout << x << " " << y << " " << z << endl;
    }
    else if(flag == 6)
    {
        int x = x_ + c->GetHigh();
        int y = y_ + c->GetWidth();
        int z = z_ + c->GetLength();
        Point cargop2(x, y, z);
        c->SetPoint2(cargop2);
        //cout << x << " " << y << " " << z << endl;
    }

}

json calcGoodList(vector<Cargo*> vec, int& index) {
	json cargoes = json::array();

	for (Cargo* cargo : vec) {
		Point bottom_left = cargo->GetPoint1();
		cargoes.push_back({
				{"materialCode", cargo->materialcode},
				{"restrictionFlag", cargo->flag},
				{"x", bottom_left.GetPointX() * 10 },
				{"y", bottom_left.GetPointY() * 10},
				{"z", bottom_left.GetPointZ() * 10},
				{"trainIndex", index++}
			});
	}

	return cargoes;
}

json calcStepList(const map<string, vector<Cargo*>>& orders, int& index)
{
	int step_index = 1;
	json steps = json::array();;

	for (auto&& pair : orders) {
		steps.push_back({
				{"step", step_index++},
				{"qty", pair.second.size()},
				{"directionNum", "1*1*1"},
				{"orderCode", pair.first},
				{"goodList", calcGoodList(pair.second, index)}
			});
	}

	return steps;
}

int calcGoodNum(const map<string, vector<Cargo*>>& orders) {
	int num = 0;

	for (auto&& pair : orders) {
		num += pair.second.size();
	}

	return num;
}

json calcTrains(const map<int, map<string, vector<Cargo*>>>& containers) {
	json trainList = json::array();

	for (auto&& pair : containers) {
		int index =0;

		int64_t totalCapacityCM3 = 0;
		double totalWeight = 0;

		for (auto&& p2 : pair.second) {
			for (Cargo* cargo : p2.second) {
				totalCapacityCM3 += cargo->GetVolumeCM3();
				totalWeight += cargo->GetWeight();
			}
		}

		double totalCapacityM3 = (double)totalCapacityCM3 / 10e6;

		Bin* currentBin = bin[pair.first];

		double packingRate = (double)totalCapacityCM3 / currentBin->GetVolumeCM3();

		trainList.push_back({
				{ "train", pair.first},
				{"modelCode", bin[pair.first]->modelCode },
				{"goodNum", calcGoodNum(pair.second)},
				{"totalCapacity", totalCapacityM3 },  // !!!
				{"totalWeight", totalWeight}, /// !!!
				{"packingRate", packingRate}, /// !!!s
				{"stepList", calcStepList(pair.second, index)}
			});
	}

	return trainList;
}

string another_ANS(vector<Cargo*>  vec, bool valid)
{
	sort(vec.begin(), vec.end(), ComparePacking2);

	json reply;

	reply["tradeId"] = global_tradeId;

	if (!valid) {
		reply["status"] = "0";
		return reply.dump();
	}

	map<int, map<string, vector<Cargo*>>> containers;

	for (Cargo* cargo : vec)
		containers[cargo->numofpack][cargo->orderCode].push_back(cargo);

	reply["status"] = "1";
	reply["carNum"] = containers.size();
	reply["goodNum"] = vec.size();

	{
		int64_t totalCapacityCM3 = 0;
		double totalWeight = 0;

		for (Cargo* cargo : vec) {
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

// string ANS(vector<Cargo*> vec, bool valid)
// {
// 	sort(vec.begin(), vec.end(), ComparePacking2);

//     int step_ = 1;
//     int numofpack = -1;

//     json all;
// 	all["tradeId"] = global_tradeId;

// 	if (valid) {
// 		all["status"] = "ok";
// 	} else {
// 		all["status"] = "failed";
// 		return all.dump();
// 	}

// 	all["msg"] = "";

//     all["carNum"] = maxnbin;
//     all["goodNum"] = maxn;
//     all["totalCapacity"] = 0;
//     all["totalWeight"] = 1;

//     json em = json::array();

//     json train[maxnbin_];
//     json trainem = json::array();
//     trainem.clear();

//     int k = 0;
//     for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
//     {
//         Cargo* c = *it;
//         if(c->numofpack != numofpack)
//         {
//             if(k != 0)
//             {
//                 train[numofpack]["stepList"] = trainem;
//                 em.push_back(train[numofpack]);
//                 trainem.clear();
//             }
//             k++;
// 			numofpack = c->numofpack;

//             train[c->numofpack]["train"] = c->numofpack;
//             train[c->numofpack]["modelCode"] = bin[c->numofpack]->modelCode;
//             train[c->numofpack]["goodNum"] = 0;
//             train[c->numofpack]["totalCapacity"] = 0;
//             train[c->numofpack]["totalWeight"] = 0;
//             train[c->numofpack]["packingRate"] = 0;
//         }

//         json st;
//         st["step"] = step_;
//         step_++;
//         st["qty"] = 1;
//         st["directionNum"] = "1*1*1";

//         json good;
//         good["materialCode"] = c->materialcode;
//         good["restrictionFlag"] = c->flag;
//         good["x"] = (c->GetPoint1()).GetPointX() * 10;
//         good["y"] = (c->GetPoint1()).GetPointY() * 10;
//         good["z"] = (c->GetPoint1()).GetPointZ()* 10;
//         good["orderCode"] = c->orderCode;
//         good["trainIndex"] = c->numofpack;

//         st["goodList"] = good;
//         trainem.push_back(st);
//     }

//     train[numofpack]["stepList"] = trainem;
//     em.push_back(train[numofpack]);

//     all["trainList"] = em;

//     //cout << all << endl;
//     string ans = all.dump();
//     return ans;
//     /*
//     int step_ = 1;
//     int numofpackL = -1;
//     int numofpackR = -1;
//     string str = "{\"status\":1,";
//     str += "\"msq\":\"\",";
//     str += "\"carNum\":";
//     str = str + to_string(maxnbin) + ",";
//     str = str + "\"goodNum\":" + to_string(maxn) + ",";
//     str = str + "\"totalCapacity\": 0,";
//     str = str + "\"totalWeight\":1,";
//     str = str + "\"trainList\":[";

//     for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
//     {
//         Cargo* c = *it;

//         if(c->numofpack != numofpackL)
//         {
//             if(it != vec.begin())
//             {
//                 str += "]}";
//             }
//             str += "{";
//             str += "\"train\":" + to_string(c->numofpack) +",";
//             str += "\"modelCode\":\"" + bin[c->numofpack]->modelCode + "\",";
//             str += "\"goodNum\":0,";
//             str += "\"totalCapacity\":0,";
//             str += "\"totalWeight\":0,";
//             str += "\"packingRate\":0,";
//             str += "\"stepList\":[";
//             numofpackL = c->numofpack;
//         }

//         str += "{\"step\":" + to_string(step_) + ",";
//         step_++;
//         str += "\"qty\":1,";
//         str += "\"directionNum\":\"1*1*1\",";
//         str += "\"goodList\":[{";

//         str += "\"materialCode\":\"" + c->materialcode + "\",";
//         str += "\"restrictionFlag\":" + to_string(c->flag) + ",";
//         str += "\"x\":" + to_string((c->GetPoint1()).GetPointX()) + ",";
//         str += "\"y\":" + to_string((c->GetPoint1()).GetPointY()) + ",";
//         str += "\"z\":" + to_string((c->GetPoint1()).GetPointZ()) + ",";
//         str += "\"orderCode\":" + to_string(c->orderCode) + ",";
//         str += "\"trainIndex\":" + to_string(c->numofpack);
//         str += "}]}";

//         if(it == vec.end() - 1)
//         {
//             str += "]}";
//         }
//     }

//     str += "]}";
//     return str;
//     */
// }

void ShowCargo(vector<Cargo*> vec)     //装箱箱子位置输出
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
    cout << str;
    const char *ch = str.c_str();
    SocketUDP(ch);

    //p2_
    for(vector<Cargo*>::iterator it = vec.begin(); it != vec.end(); it++)
    {
        Cargo* c = *it;
        printf("输出货物的位置:\n");
        ShowPoint(c->GetPoint1(), c->GetPoint2());
        cout << "货物编号：" << c->GetId() << " 订单号: " << c->ding << " 被装集装箱序列号：" << c->numofpack << endl;
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
        cout << str;
        const char *ch = str.c_str();
        if((*it)->numofpack == CHE && ONE == 1)
        {
            //Sleep(1000);
            SocketUDP(ch);
        }
    }
}
//更新最大空闲空间
void StateUpdate(Cargo *c, int flag)
{
    //每次装箱子将之前的空间给切碎并重新维护
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
        //证明两个立方体相交
        if(x1 < x2 && y1 < y2 && z1 < z2)
        {
            Point p1(x1, y1, z1), p2(x2, y2, z2);
            //大的立方体是ems(X1-X2)_, 小的立方体是p1-p2(X3-X4)
            int X1 = (ems_.GetPoint1()).GetPointX(), Y1 = (ems_.GetPoint1()).GetPointY(), Z1 = (ems_.GetPoint1()).GetPointZ();
            int X2 = (ems_.GetPoint2()).GetPointX(), Y2 = (ems_.GetPoint2()).GetPointY(), Z2 = (ems_.GetPoint2()).GetPointZ();
            int X3 = x1, Y3 = y1, Z3 = z1;
            int X4 = x2, Y4 = y2, Z4 = z2;
            //切为6个空间
            //1(左侧)
            if(X1 != X3 && Y1 != Y2 && Z1 != Z2)
            {
                Point point1(X1, Y1, Z1), point2(X3, Y2, Z2);
                int length[3];
                length[0] = X3 - X1;    //最小的长度
                length[1] = Y2 - Y1;
                length[2] = Z2 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;
                    vece2.push_back(ems);
                }
            }
            //2(右侧)
            if(X4 != X2 && Y1 != Y2 && Z1 != Z2)
            {
                Point point1(X4, Y1, Z1), point2(X2, Y2, Z2);
                int length[3];
                length[0] = X2 - X4;    //最小的长度
                length[1] = Y2 - Y1;
                length[2] = Z2 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;
                    vece2.push_back(ems);
                }
            }
            //3（前侧）
            if(X1 != X2 && Y1 != Y3 && Z1 != Z2)
            {
                Point point1(X1, Y1, Z1), point2(X2, Y3, Z2);
                int length[3];
                length[0] = X2 - X1;    //最小的长度
                length[1] = Y3 - Y1;
                length[2] = Z2 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;
                    vece2.push_back(ems);
                }
            }
            //4（后侧）
            if(X1 != X2 && Y4 != Y2 && Z1 != Z2)
            {
                Point point1(X1, Y4, Z1), point2(X2, Y2, Z2);
                int length[3];
                length[0] = X2 - X1;    //最小的长度
                length[1] = Y2 - Y4;
                length[2] = Z2 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;
                    vece2.push_back(ems);
                }
            }
            //5（下侧）
            if(X1 != X2 && Y1 != Y2 && Z1 != Z3)
            {
                Point point1(X1, Y1, Z1), point2(X2, Y2, Z3);
                int length[3];
                length[0] = X2 - X1;    //最小的长度
                length[1] = Y2 - Y1;
                length[2] = Z3 - Z1;
                sort(length, length+3);
                if(length[0] > minCargo1 && length[1] > minCargo2 && length[2] >= minCargo3)
                {
                    EMSpace ems(point1, point2);
                    ems.level = (*ite).level;
                    ems.yes = (*ite).yes;
                    ems.weight = (*ite).weight;
                    vece2.push_back(ems);
                }
            }
            //6（上侧）
            if(X1 != X2 && Y1 != Y2 && Z4 != Z2)
            {
                Point point1(X1, Y1, Z4), point2(X2, Y2, Z2);
                int length[3];
                length[0] = X2 - X1;    //最小的长度
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
                    vece2.push_back(ems);
                }
            }
        }
        //两个立方体不相交
        else
        {
            vece2.push_back(ems_);
        }
    }
    //vece中存储为存在冗余的空间
    vece.clear();
    for(ite2 = vece2.begin(); ite2 != vece2.end(); ite2++)
    {
        vece.push_back(*ite2);
    }
    //将vece中存储的冗余空间删除掉
    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        //ShowPoint(ite->GetPoint1(), ite->GetPoint2());
        for(ite2 = vece2.begin(); ite2 != vece2.end(); ite2++)
        {
            if(Compare(ite->GetPoint1(), ite2->GetPoint1()) == 2 && Compare(ite->GetPoint2(), ite2->GetPoint2()) == 2)
            {
                continue;
            }
            //ite包含在ite2中
            else if((Compare(ite->GetPoint1(), ite2->GetPoint1()) == 1 || Compare(ite->GetPoint1(), ite2->GetPoint1()) == 2) && (Compare(ite->GetPoint2(), ite2->GetPoint2()) == -1 || Compare(ite->GetPoint2(), ite2->GetPoint2()) == 2))
            {
                vece.erase(ite);
                ite--;
                break;
            }
        }
    }

    //将vece中里面的空间删除掉
    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        Point pv_1 = ite->GetPoint1();  //空间的两点
        Point pv_2 = ite->GetPoint2();
        Point pc_1 = c->GetPoint1();    //箱子的两点
        Point pc_2 = c->GetPoint2();
        if(pv_2.GetPointX() <= pc_1.GetPointX() && pv_1.GetPointZ() <= pc_2.GetPointZ())
        {
            vece.erase(ite);
            ite--;
        }
    }

    for(ite = vece.begin(); ite != vece.end(); ite++)
    {
        Point pv_1 = ite->GetPoint1();  //空间的两点
        Point pv_2 = ite->GetPoint2();
        Point pc_1 = c->GetPoint1();    //箱子的两点
        Point pc_2 = c->GetPoint2();

        if(pv_1.GetPointX() + LengthOfPacking <= pc_1.GetPointX())      //避免箱子放置过于深入
        {
            //cout << pv_1.GetPointX() << "    " <<  pc_1.GetPointX() << "    ";
            int x_ = pc_2.GetPointX();
            //cout << x_ - LengthOfPacking << "   ";
            ite->Change(x_ - LengthOfPacking);
            //cout << (ite->GetPoint1()).GetPointX() << endl;
        }

    }
}

template <typename T>
T field_of(const json& object, const char* name) {
	try {
		return T(object.at(name));
	} catch (const exception& err) {
		throw runtime_error(string("getting `") + name + string("`: ") + err.what());
	}
}

void JSON(string read)
{

    json input = json::parse(read);//输入json对象


    // string tradeId = field_of<string>(input, "tradeId"); //tradeId号
	global_tradeId = field_of<string>(input, "tradeId"); //tradeId号

    json orderList = field_of<json>(input, "orderList");//订单列表
    json binList = field_of<json>(input, "vehicleModelList");//货箱列表
    int cargoid = 0;

    //cout << orderList.size() << endl;
    for(int i=0; i<orderList.size(); i++) //遍历所有订单
    {
        json order = orderList[i];//当前订单

        // long long orderCode = order.at("orderCode");//订单ID,string用来中转
        string orderCode = field_of<string>(order, "orderCode");//订单ID,string用来中转

        //int orderID = stoi(tempString);
        //cout << orderID << endl;
        int unloadingSequence = field_of<int>(order, "unloadingSequence");

        json boxList = field_of<json>(order,"goodList");//货物集合json对象
        //cout << boxList.size() << endl;


        for(int j =0; j<boxList.size(); j++) //遍历当前订单的所有货物
        {

            json box = boxList[j];//当前货物

            //string qty = box.at("qty"); //货物个数
            //int qtynum = stoi(qty);
            string materialCode = field_of<string>(box,"materialCode");//订单ID,string用来中转
            int qtynum = field_of<int>(box, "qty");
            //cout << qtynum << endl;
            qtynum--;

            string idOfSet = field_of<string>(box,"setCode");//套机编码
            //cout << idOfSet << endl;


            //string tempLength;//得到长宽高和重量，string用来中转
            /*
            tempLength=box.at("length");
            int x = stoi(tempLength);
            cout << x << endl;
            tempLength=box.at("width");
            int y = stoi(tempLength);
            cout << y << endl;
            tempLength=box.at("height");
            int z = stoi(tempLength);
            cout << z << endl;
            tempLength=box.at("weight");
            int weight = stoi(tempLength);
            cout << weight << endl;
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
            ca[cargoid]->ding = unloadingSequence;
            ca[cargoid]->weight = weight;
            ca[cargoid]->materialcode = materialCode;
            ca[cargoid]->orderCode = orderCode;
            //cout << ca[cargoid]->weight << endl;

            json boxLimit = field_of<json>(box, "restrictionList");//装箱限制json对象
            for(int k =0; k<boxLimit.size(); k++) //遍历所有限制
            {
                json oriLimit = boxLimit[k];//当前限制规则

                int orientation = stoi(field_of<string>(oriLimit, "flag"));//限制方向
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
                //cout << orientation << endl;
                //cout << "     " << direction << endl;

                bool tempIsBearSurface = field_of<bool>(oriLimit,"isBear");//承重面限制
                int isBearSurface=0;
                int bearLevel = 1;

                if(tempIsBearSurface==true)
                {
                    isBearSurface=1;
					bearLevel = field_of<int>(oriLimit, "bearLevel");
                }

                //cout << isBearSurface << endl;

                //string tempBearLevel = oriLimit.at("bearLevel");//承重级别限制
                //int bearLevel = stoi(tempBearLevel);

                //cout << bearLevel << endl << endl;

                //bool isStackLimit = oriLimit.at("isStack");//堆码限制
                //int stackLevel = 0;//堆码层数
                //int stackLevel = oriLimit.at("stackLevel");

                ca[cargoid]->SetOri(direction, 0);
                (ca[cargoid]->node[direction]).l = bearLevel;
                (ca[cargoid]->node[direction]).y = isBearSurface;
                (ca[cargoid]->node[direction]).f = 1;
            }
            vecc.push_back(ca[cargoid]);
            cargoid++;

            while(qtynum)
            {
                //cout << qtynum << endl;
                ca[cargoid] = new Cargo(x, y, z);
                ca[cargoid]->Setnum(0);
                for(int j=0; j<6; j++)
                {
                    ca[cargoid]->SetOri(j, 0);
                    (ca[cargoid]->node[j]).l = (ca[cargoid-1]->node[j]).l;
                    (ca[cargoid]->node[j]).f = (ca[cargoid-1]->node[j]).f;
                    (ca[cargoid]->node[j]).y = (ca[cargoid-1]->node[j]).y;
                }
                ca[cargoid]->SetId(cargoid);
                ca[cargoid]->ding = unloadingSequence;
                ca[cargoid]->weight = weight;
				ca[cargoid]->orderCode = orderCode;

				ca[cargoid]->materialcode = materialCode;

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

    for(int i=0; i<binList.size(); i++) //遍历货箱列表
    {

        json binq = binList[i];//当前货箱

        //json binAttribute = bin.at("货箱属性");//货箱属性

        //string idOfBin = binAttribute.at("货箱型号");//货箱型号
        string typeOfBin = field_of<string>(binq,"modelCode");//货箱类型

        /*
        string intTemp=bin.at("weight");//转int临时变量
        int weightOfBin = stoi(intTemp);//货箱重量

        intTemp = bin.at("length");
        int x = stoi(intTemp);//货箱长

        intTemp = bin.at("width");
        int y = stoi(intTemp);//货箱宽

        intTemp = bin.at("height");
        int z = stoi(intTemp);//货箱高

        string binNumbers = bin.at("qty");//货箱可用数量
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

    maxn = cargoid;
    maxnbin = binid;
}



char* binPackingAlgorithm(const char* str)
{
    Init();
    JSON(str);

    Chromosome* chromosome_ans = GeneticAlgorithm();
    bool valid = chromosome_ans && chromosome_ans->PackingOf(vecc);

    string ans = another_ANS(vecc, valid);

	char* reply = (char*)malloc(ans.size());
	memcpy(reply, ans.data(), ans.size());
 // strcpy(reply, ans.c_str());

	delete_vector(bin, maxnbin);

	delete_vector(ca, maxn);

	for (auto ptr : allocated_chrosomes)
		delete ptr;

	return reply;
}

extern "C" JNIEXPORT jstring JNICALL Java_com_haier_tms_jni_NativeAlgorithm_binPackingAlgorithm
  (JNIEnv * env, jclass cls, jstring str)
{
	const char *nativeString = env->GetStringUTFChars(str, 0);
	char* reply =	binPackingAlgorithm(nativeString);

	env->ReleaseStringUTFChars(str, nativeString);
	jstring reply_jstring = env->NewStringUTF(reply);
	free(reply);
	return reply_jstring;
}

extern "C" JNIEXPORT void JNICALL Java_com_haier_tms_jni_NativeAlgorithm_hello
(JNIEnv *, jclass)
{
	ofstream os("/root/hello_output.txt", std::ios_base::out);
    os << "hello world from c world" << endl;

}


int main()
{
    Init();

	ostringstream os;
	os << cin.rdbuf();

	char* str = binPackingAlgorithm(os.str().c_str());
	cerr << str << endl;
	free(str);
}
