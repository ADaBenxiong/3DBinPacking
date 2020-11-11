#ifndef POINT_H_INCLUDED
#define POINT_H_INCLUDED

// #include<bits/stdc++.h>
#include <vector>
using namespace std;
class Point
{
public:
    Point(){}
    Point(int x_, int y_, int z_)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    ~Point(){}
    int GetPointX() const;
    int GetPointY() const;
    int GetPointZ() const;
private:
    int x;
    int y;
    int z;
};

int Point::GetPointX() const
{
    return x;
}

int Point::GetPointY() const
{
    return y;
}

int Point::GetPointZ() const
{
    return z;
}

#endif // POINT_H_INCLUDED
