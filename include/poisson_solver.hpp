#include <vector>
#include <string>

int idx(int i, int j, int nx);
// 输入一个矩阵，矩阵中一个元素的位置，得到存储在vector中的位置

void initial_t_field(std::vector<double> &u,
                     double t_top,
                     double t_bottom,
                     double t_left,
                     double t_right,
                     int nx,
                     int ny);
// 初始化一个温度场，需要输入一个初始化为0的vector，以及四个边界温度，直接引用目标vector

void write_field(const std::vector<double> &u, const std::string &filename, int nx, int ny);
// 输入一个温度场，文件名，以及网格长度宽度，将温度场写入到文件名中

void update_field_single(std::vector<double> &u_old,
                         std::vector<double> &u_new,
                         int nx,
                         int ny);
// 用来按照中心点温度是四周的平均值方式更新温度,仅仅更新一次

void update_field_final(std::vector<double> &u_old,
                        std::vector<double> &u_new,
                        int nx,
                        int ny,
                        double stop_dif);
// 用来按照中心点温度是四周的平均值方式更新温度,直至更新前后对应点的最大温差绝对值小于 stop_dif
