#include <fstream>
#include <iostream>

#include "json.hpp"

using namespace std;
using json = nlohmann::json;

string JSON()
{
    ifstream read("input structure.json");

	json input = json::parse(read);//输入json对象


	json orderList = input.at("订单列表");//订单列表
	json binList = input.at("货箱列表");//货箱列表

	for(int i=0;i<orderList.size();i++){//遍历所有订单

		json order = orderList[i];//当前订单

		string tempString = order.at("订单ID");//订单ID,string用来中转
		int orderID = stoi(tempString);
		cout << orderID << endl;

		json boxList = order.at("货物集合");//货物集合json对象

		if(boxList.size()==1){//当货物集合大小为1时即当前没有货物(为了解决示例中的省略号情况，实际应该是需要删了)
			continue;
		}

		for(int j =0;j<boxList.size();j++){//遍历当前订单的所有货物

			json box = boxList[j];//当前货物

			string idOfSet = box.at("套机编码");//套机编码
			cout << idOfSet << endl;

			string tempLength;//得到长宽高和重量，string用来中转
			tempLength=box.at("长");
			int x = stoi(tempLength);
			cout << x << endl;
			tempLength=box.at("宽");
			int y = stoi(tempLength);
			cout << y << endl;
			tempLength=box.at("高");
			int z = stoi(tempLength);
			cout << z << endl;
			tempLength=box.at("重量");
			int weight = stoi(tempLength);
			cout << weight << endl;

			json boxLimit = box.at("装箱限制");//装箱限制json对象
			for(int k =0;k<boxLimit.size();k++){//遍历所有限制

				json oriLimit = boxLimit[k];//当前限制规则

				string orientation = oriLimit.at("摆放方向");//限制方向
				int direction = 0;
				if(orientation == "立放正向")
                    direction = 1;
                else if(orientation == "立放横向")
                    direction = 3;
                else if(orientation == "侧放正向")
                    direction = 2;
                else if(orientation == "侧放横向")
                    direction = 4;
                else if(orientation == "卧放正向")
                    direction = 5;
                else if(orientation == "卧放横向")
                    direction = 6;
                //cout << orientation << endl;
                cout << "     " << direction << endl;

				string tempIsBearSurface = oriLimit.at("是否承重面");//承重面限制
				bool isBearSurface=false;
				if(tempIsBearSurface=="是"){
					isBearSurface=true;
				}

				cout << isBearSurface << endl;

				string tempBearLevel = oriLimit.at("承重级别");//承重级别限制
				int bearLevel = stoi(tempBearLevel);

				cout << bearLevel << endl << endl;

				string isStackLimit = oriLimit.at("堆码限制");//堆码限制
				int stackLevel = 0;//堆码层数
				if(isStackLimit=="是"){
					string tempStackLevel = oriLimit.at("堆码层数");
					stackLevel = stoi(tempStackLevel);
				}
			}
		}
	}

	for(int i=0;i<binList.size();i++){//遍历货箱列表

		json bin = binList[i];//当前货箱

		json binAttribute = bin.at("货箱属性");//货箱属性

		string idOfBin = binAttribute.at("货箱型号");//货箱型号
		string typeOfBin = binAttribute.at("货箱类型");//货箱类型

		string intTemp=binAttribute.at("重量");//转int临时变量
		int weightOfBin = stoi(intTemp);//货箱重量

		intTemp = binAttribute.at("长");
		int x = stoi(intTemp);//货箱长

		intTemp = binAttribute.at("宽");
		int y = stoi(intTemp);//货箱宽

		intTemp = binAttribute.at("高");
		int z = stoi(intTemp);//货箱高

		string binNumbers = bin.at("可用数量");//货箱可用数量
		int numOfBin = stoi(binNumbers);

	}

	read.close();
}
