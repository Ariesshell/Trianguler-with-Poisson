#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "poissonDisc.h"
using namespace std;
#define PI       3.14159265358979323846

bool inCircle(P a, P b, double r) {
    a.first -= b.first;
    a.second -= b.second;
    return a.first * a.first + a.second * a.second < r* r;
}

bool isValid(double x, double y, double r, double cellSize, P newP, vector<vector<int>>& grid, vector<P>& points) {
    //if (newP.first<0 || newP.first>h || newP.second<0 || newP.second>w)
      //  return false; //出界

    // 计算出有可能发生碰撞五乘五格子的上下左右边界
    int up = max(0, int(newP.first / cellSize) - 2);
    int down = min(int(grid.size() - 1), int(newP.first / cellSize) + 2);

    int left = max(0, int(newP.second / cellSize) - 2);
    int right = min(int(grid[0].size() - 1), int(newP.second / cellSize) + 2);

    for (int i = up; i <= down; ++i) {
        for (int j = left; j <= right; ++j) {
            if (grid[i][j] < 0) {
                continue;
            }
            if (inCircle(newP, points[grid[i][j]], r)){
                return false;
            }
        }
    }
    return true;
}

bool isOpacity(double w, double h, P newP, unsigned char* data) {
    int x = (int)newP.first;
    int y = (int)newP.second;
    y = max(0, y - 1);
    int alpha = data[(y * (int)w + x)*4]+ data[(y * (int)w + x+1) * 4] + data[(y * (int)w + x+2) * 4];
    //printf("index:%d\n", data[y * (int)w * 4 + x] + data[y * (int)w * 4 + x + 1] + data[y * (int)w * 4 + x + 2]);
    return alpha > 0 ? true : false;
}

vector<P> sampling(double h, double w, double r, unsigned char* data, int numSamplesBeforeReject) {
    srand(NULL);
    double cellSize = r / sqrt(2);
    vector<P> samplePoints;
    vector<P> processPoints;
    // 划分存储网格
    vector<vector<int>> grid(ceil(h / cellSize), vector<int>(ceil(w / cellSize), -1));
    // 初始种子点
    processPoints.push_back({ h / 2, w / 2 });
    processPoints.push_back({ 0, 0 });
    processPoints.push_back({ 0, w });
    processPoints.push_back({ h ,0 });
    processPoints.push_back({ h ,w });

    while (!processPoints.empty()) {
        bool candidateAccepted = false;
        // 从processPoints中随机取一个点
        int processIdx = rand() % processPoints.size();
        auto newP = processPoints[processIdx];

        // 每次按种子为中心，从半径r到2r的地带随机选取一点。
        // 并判断是否可行，可行则退出，并加入到候选生成种子，否继续，直到采样次数用完。
        for (int i = 0; i < numSamplesBeforeReject; ++i) {
            double angle = rand() % 360 * 1.0 / 360 * 2 * PI;
            auto dir = P(sin(angle), cos(angle));
            double rl = rand() % 1000 * 1.0 / 1000 * r;
            newP.first += dir.first * (rl + r);
            newP.second += dir.second * (rl + r);

            //if (newp.first<0 || newp.first>h || newp.second<0 || newp.second>w)
            //    continue; //出界
            newP.first = max((float)0, newP.first);
            newP.first = min((float)h, newP.first);
            newP.second = max((float)0, newP.second);
            newP.second = min((float)w, newP.second);

            bool res = isValid(h, w, r, cellSize, newP, grid, samplePoints);
            res &= isOpacity(h, w, newP, data);
            if (res) {
                candidateAccepted = true;
                grid[newP.first / cellSize][newP.second / cellSize] = samplePoints.size();
                samplePoints.push_back(newP);
                processPoints.push_back(newP);
                break;
            }
        }

        // 采样次数用完，没有得到可行结果，说明该种不可用，直接删除。
        if (!candidateAccepted) {
            if (processIdx != processPoints.size() - 1)
                processPoints[processIdx] = processPoints[processPoints.size() - 1];

            processPoints.pop_back();
        }
    }

    return samplePoints;
}