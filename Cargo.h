#ifndef CARGO_H_INCLUDED
#define CARGO_H_INCLUDED

#include "Point.h"
#include <string>
#include <cstdint>
#include <limits>

//表示用于装箱的货物
class Cargo
{
public:
    Cargo(){}
    Cargo(int x_, int y_, int z_)
    {
        x = x_;
        y = y_;
        z = z_;
        for(int i =0; i<6; i++)
        {
            node[i].f = 0;
            node[i].l = 0;
            // if(i == 0 || i == 2)
            //     node[i].y = 1;
            // else
			node[i].y = 0;
			node[i].isStack = false;
			node[i].stackLevel = 100;
        }
        numofpack = -1;
        yesk = 0;
    }
    ~Cargo(){}
    //确定箱子的两点坐标
    void SetPoint1(Point p)
    {
        p1 = p;
    }
    void SetPoint2(Point p)
    {
        p2 = p;
    }
    Point GetPoint1() const
    {
        return p1;
    }
    Point GetPoint2() const
    {
        return p2;
    }
    //得到箱子的长宽高
    int GetLength() const
    {
        return x;
    }
    int GetWidth() const
    {
        return y;
    }
    int GetHigh() const
    {
        return z;
    }

	// cm^3
	int64_t GetVolumeCM3() const {
		return x * y * z;
	}

	double GetWeight() const {
		return weight;
	}

    //箱子放置方向选择
    int Getnum() const
    {
        return num;
    }
    void Setnum(int k)
    {
        num = k;
    }
    int GetOri(int i) const
    {
        return ori[i];
    }
    void SetOri(int i, int k)
    {
        ori[i] = k;
    }
    //解码染色体后获得bps，vbo, BPS是序列, VBO是方向
    double GetVBO() const
    {
        return vbo;
    }
    double GetBPS() const
    {
        return bps;
    }
    void SetVBO(double k)
    {
        vbo = k;
    }
    void SetBPS(double k)
    {
        bps = k;
        if(yesk == 1)
        {
            //c->SetBPS_(k + 0.001);
        }
    }
    void SetBPS_(double k)
    {
        bps = k;
    }
    int GetId() const
    {
        return id;
    }
    void SetId(int id_)
    {
        id = id_;
    }
    struct Node
    {
        int f;//方向
        int y;//称重面
        int l;//承重级别
		bool isStack;	//堆码限制
		int stackLevel;	//堆码层数
    };
    Node node[6];
    int ding;   //订单序列号
    int numofpack;  //被装集装箱号
    double weight; //表示箱子的重量

    int yesk;   //表示是否是套机
    char str[10];   //套机编码
    Cargo* c;   //套机

    int flag = 1;

	string materialCode;
    // long long orderCode;

    string orderCode;
	string set_code;

private:
    //箱子的长宽高
    int x;
    int y;
    int z;

    int id;  //表示箱子的编号
    double bps = 0; //染色体解码后第一个值
    double vbo = 0;//染色体解码后第二个值

    //箱子被装箱后的两点坐标
    Point p1;
    Point p2;
    int num;    //可以放置方向个数
    int ori[6]; //方向向量（分别为长宽高）  x,y,z设置方向1， x,z,y设置方向2， y,x,z设置方向3, y,z,x设置方向4, z,x,y设置方向5, z,y,x设置方向6
};


class EMSpace
{
public:
    EMSpace(Point x, Point y)
    {
        p1 = x;
        p2 = y;
        yes = 1;
        level = 0;      //承重级别为level
        weight = 300;

		layer = 100;
		max_weight = numeric_limits<double>::max();
    }
    ~EMSpace(){}
    //得到两个顶点坐标
    Point GetPoint1() const
    {
        return p1;
    }
    Point GetPoint2() const
    {
        return p2;
    }
    //得到长宽高
    int GetLength() const
    {
        return p2.GetPointX() - p1.GetPointX();
    }
    int GetWidth() const
    {
        return p2.GetPointY() - p1.GetPointY();
    }
    int GetHigh() const
    {
        return p2.GetPointZ() - p1.GetPointZ();
    }
    void Change(int p)
    {
        Point pnew(p, p1.GetPointY(), p1.GetPointZ());
        p1 = pnew;
    }
    int level;  //表示这个空间可以放的称重级别
    double weight; //表示这个空间箱子的重量
    int yes;    //是否是称重面

	int layer;	//表示该空间最多堆码几层
	double max_weight;	//表示该空间最大承重
	string cargo_code;	//表示空间下面堆放货物的编码

private:
    Point p1;
    Point p2;
};

class Bin
{
public:
    Bin(Point x, Point y)
    {
        p1 = x;
        p2 = y;
    }
    ~Bin(){}
    void SetPoint(Point p1_, Point p2_)
    {
        p1 = p1_;
        p2 = p2_;
    }
    //两个顶点坐标
    Point GetPoint1() const
    {
        return p1;
    }
    Point GetPoint2() const
    {
        return p2;
    }
    //箱子编号的问题
    void SetNum(int k)
    {
        num = k;
    }
    int GetNum() const
    {
        return num;
    }
    //得到长宽高
    int64_t GetLength() const
    {
        return p2.GetPointX() - p1.GetPointX();
    }
    int64_t GetWidth() const
    {
        return p2.GetPointY() - p1.GetPointY();
    }
    int64_t GetHigh() const
    {
        return p2.GetPointZ() - p1.GetPointZ();
    }

	int64_t GetVolumeCM3() const {
		return GetLength() * GetHigh() * GetWidth();
	}

    string modelCode;
private:
    Point p1;
    Point p2;
    int num;    //表示箱子的编号（对于多箱的问题）
};


#endif // CARGO_H_INCLUDED
