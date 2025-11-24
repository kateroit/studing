#define GLFW_INCLUDE_NONE
#include <iostream>
#include <filesystem>
#include <fstream>
#include <math.h>
#include <vector>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

//глобальные переменные для управления камерой
GLfloat d[3] = { 0.0, 0.0, -0.8 };
GLfloat viewangle = 0;
GLfloat tippangle = 0;
GLfloat scaleF = 0.2;

//оси координат
float XUP[3] = { 1,0,0 }, XUN[3] = { -1, 0, 0 },
YUP[3] = { 0,1,0 }, YUN[3] = { 0,-1, 0 },
ZUP[3] = { 0,0,1 }, ZUN[3] = { 0, 0,-1 },
ORG[3] = { 0,0,0 };

//прототипы функций
static void keys(GLFWwindow* window, int key, int scancode, int action, int mods);
void normalize(float* v);
void calculate_vector(float* a, float* b, float* result);
void calculate_normal(float* a, float* b, float* c, float* result);
void writeVRML(const std::vector<double>& vertices_x, const std::vector<double>& vertices_y, const std::vector<double>& vertices_z, const std::vector<int>& faces);
GLFWwindow* initGLFW();
bool loadDepthMap(const std::string& filePath, std::vector<std::vector<double>>& dmap, int& height, int& width);
void convertTo3DCoordinates(const std::vector<std::vector<double>>& dmap, std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz);
void createTriangleMesh(int height, int width, const std::vector<std::vector<double>>& z, std::vector<int>& faces);
void drawCoordinateAxes();
void drawModel(const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz, const std::vector<int>& faces);
void applyTransformations();
void renderFrame(GLFWwindow* window, const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz, const std::vector<int>& faces);
void printControls();
void exportToVRML(const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz, const std::vector<int>& faces);
void cleanup(GLFWwindow* window);

//нормализация вектора 
void normalize(float* v)
{
    float length = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    for (int i = 0; i < 3; i++)
    {
        v[i] = v[i] / length;
    }
}

//вычисление векторного произведения
void calculate_vector(float* a, float* b, float* result)
{
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = -(a[0] * b[2] - a[2] * b[0]);
    result[2] = a[0] * b[1] - a[1] * b[0];
    normalize(result);
}

//рассчет нормали 
void calculate_normal(float* a, float* b, float* c, float* result)
{
    float x[] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
    float y[] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
    calculate_vector(x, y, result);
}

//перемещение, вращение и тд засчет клавиш 
static void keys(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS) {
        switch (key) {
        case GLFW_KEY_ESCAPE:
            glfwSetWindowShouldClose(window, GL_TRUE);
            break;
        case GLFW_KEY_Z:
            scaleF += 0.1;
            break;
        case GLFW_KEY_X:
            scaleF -= 0.1;
            break;
        case GLFW_KEY_W:
            tippangle -= 5;
            break;
        case GLFW_KEY_S:
            tippangle += 5;
            break;
        case GLFW_KEY_A:
            viewangle -= 5;
            break;
        case GLFW_KEY_D:
            viewangle += 5;
            break;
        case GLFW_KEY_J:
            d[0] += 0.1;
            break;
        case GLFW_KEY_K:
            d[1] += 0.1;
            break;
        case GLFW_KEY_L:
            d[2] += 0.1;
            break;
        case GLFW_KEY_R:
            d[0] -= 0.1;
            break;
        case GLFW_KEY_T:
            d[1] -= 0.1;
            break;
        case GLFW_KEY_Y:
            d[2] -= 0.1;
            break;
        }
    }
}

//запись в формат VRML
void writeVRML(const std::vector<double>& vertices_x, const std::vector<double>& vertices_y,
    const std::vector<double>& vertices_z, const std::vector<int>& faces)
{
    std::ofstream Model("model.wrl");

    //запись заголовока
    Model << "#VRML V2.0 utf8\n\n";
    Model << "Shape {\n";
    Model << "  appearance Appearance {\n";
    Model << "    material Material {\n";
    Model << "      diffuseColor 0.8 0.6 0.8\n";
    Model << "    }\n";
    Model << "  }\n";
    Model << "  geometry IndexedFaceSet {\n";
    Model << "    coord Coordinate {\n";
    Model << "      point [\n";

    //запись вершин
    for (size_t i = 0; i < vertices_x.size(); i++)
    {
        Model << "        " << vertices_x[i] << " " << vertices_y[i] << " " << vertices_z[i];
        if (i < vertices_x.size() - 1)
            Model << ",";
        Model << "\n";
    }

    Model << "      ]\n";
    Model << "    }\n";
    Model << "    coordIndex [\n";

    //запись граней
    for (size_t i = 0; i < faces.size(); i++)
    {
        if (faces[i] == -1)
        {
            Model << "      -1";
            if (i < faces.size() - 1)
                Model << ",";
            Model << "\n";
        }
        else
        {
            Model << "      " << faces[i];
            if (i < faces.size() - 1 && faces[i + 1] != -1)
                Model << ",";
        }
    }

    Model << "    ]\n";
    Model << "  }\n";
    Model << "}\n";

    Model.close();
}

//инициализация библ GLFW
GLFWwindow* initGLFW()
{
    if (!glfwInit())
        return nullptr;

    GLFWwindow* window = glfwCreateWindow(800, 800, "3D Depth Map Viewer - VRML Export", NULL, NULL);
    if (window)
    {
        glfwMakeContextCurrent(window);
        glfwSetKeyCallback(window, keys);
        glEnable(GL_DEPTH_TEST);
    }
    return window;
}

//загрузка карты глубины из файла
bool loadDepthMap(const std::string& filePath, std::vector<std::vector<double>>& dmap, int& height, int& width)
{
    std::ifstream ifs(filePath, std::ios::binary);
    if (!ifs)
    {
        std::cout << "Data file not found: " << filePath << std::endl;
        std::cout << "Current path: " << std::filesystem::current_path() << std::endl;
        return false;
    }

    std::cout << "File found: " << filePath << std::endl;

    double dheight, dwidth;
    ifs.read(reinterpret_cast<char*>(&dheight), sizeof dheight);
    ifs.read(reinterpret_cast<char*>(&dwidth), sizeof dwidth);

    height = static_cast<int>(dheight);
    width = static_cast<int>(dwidth);

    std::cout << "Height: " << height << ", Width: " << width << std::endl;

    dmap.resize(height, std::vector<double>(width));
    for (auto& row : dmap)
    {
        for (double& col : row)
            ifs.read(reinterpret_cast<char*>(&col), sizeof col);
    }

    ifs.close();
    return true;
}

//преобразование карты глубины в 3D координаты
void convertTo3DCoordinates(const std::vector<std::vector<double>>& dmap,
    std::vector<double>& vx, std::vector<double>& vy, std::vector<double>& vz)
{
    double fx = 525.0;
    double fy = 525.0;
    double cx = 319.5;
    double cy = 239.5;

    int height = dmap.size();
    int width = dmap[0].size();

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            double z = dmap[i][j] / 500.0;
            double x = (j - cx) * z / fx;
            double y = (i - cy) * z / fy;

            vx.push_back(x);
            vy.push_back(y);
            vz.push_back(z);
        }
    }
}

//создание треугольной сетки из точек
void createTriangleMesh(int height, int width, const std::vector<std::vector<double>>& z, std::vector<int>& faces)
{
    for (int i = 0; i < height - 1; i++)
    {
        for (int j = 0; j < width - 1; j++)
        {
            //пропуск точек с нулевой глубиной
            if (z[i][j] != 0 && z[i][j + 1] != 0 && z[i + 1][j] != 0 && z[i + 1][j + 1] != 0)
            {
                //первый полигон
                int idx1 = j + i * width;
                int idx2 = j + (i + 1) * width;
                int idx3 = j + 1 + i * width;

                faces.push_back(idx1);
                faces.push_back(idx2);
                faces.push_back(idx3);
                faces.push_back(-1);

                //второй полигон
                int idx4 = j + 1 + i * width;
                int idx5 = j + (i + 1) * width;
                int idx6 = j + 1 + (i + 1) * width;

                faces.push_back(idx4);
                faces.push_back(idx5);
                faces.push_back(idx6);
                faces.push_back(-1);
            }
        }
    }
}

//отрисовка осей координат
void drawCoordinateAxes()
{
    //ось X (красная)
    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    glVertex3fv(ORG); glVertex3fv(XUP);
    glEnd();

    //ось Y (зеленая)
    glColor3f(0.0, 1.0, 0.0);
    glBegin(GL_LINES);
    glVertex3fv(ORG); glVertex3fv(YUP);
    glEnd();

    //ось Z (синяя)
    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    glVertex3fv(ORG); glVertex3fv(ZUP);
    glEnd();
}

//отрисовка 3D модели
void drawModel(const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz, const std::vector<int>& faces)
{
    glBegin(GL_TRIANGLES);
    glColor3f(0.74f, 0.55f, 0.76f);

    for (size_t i = 0; i < faces.size(); i += 4)
    {
        if (i + 2 < faces.size() && faces[i] != -1 && faces[i + 1] != -1 && faces[i + 2] != -1)
        {
            int idx1 = faces[i];
            int idx2 = faces[i + 1];
            int idx3 = faces[i + 2];

            if (idx1 >= 0 && idx1 < (int)vx.size() &&
                idx2 >= 0 && idx2 < (int)vx.size() &&
                idx3 >= 0 && idx3 < (int)vx.size())
            {
                float v1[3] = { (float)vx[idx1], (float)vy[idx1], (float)vz[idx1] };
                float v2[3] = { (float)vx[idx2], (float)vy[idx2], (float)vz[idx2] };
                float v3[3] = { (float)vx[idx3], (float)vy[idx3], (float)vz[idx3] };

                float normal[3];
                calculate_normal(v1, v2, v3, normal);
                glNormal3fv(normal);

                glVertex3fv(v1);
                glVertex3fv(v2);
                glVertex3fv(v3);
            }
        }
    }
    glEnd();
}

void applyTransformations()
{
    glRotatef(tippangle, 1, 0, 0);
    glRotatef(viewangle, 0, 1, 0);
    glTranslatef(d[0], d[1], d[2]);
    glScalef(scaleF, scaleF, scaleF);
}

//рендер кадраа
void renderFrame(GLFWwindow* window, const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz, const std::vector<int>& faces)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0.2f, 0.3f, 0.3f, 1.0f);

    glLoadIdentity();
    applyTransformations();

    drawCoordinateAxes();
    drawModel(vx, vy, vz, faces);

    glfwSwapBuffers(window);
    glfwPollEvents();
}

//вывод информации об управлении
void printControls()
{
    std::cout << "Controls:\n";
    std::cout << "A, D - rotation around X axis\n";
    std::cout << "W, S - rotation around Y axis\n";
    std::cout << "Z, X - zoom in/out\n";
    std::cout << "J, R - move along X axis\n";
    std::cout << "K, T - move along Y axis\n";
    std::cout << "L, Y - move along Z axis\n";
    std::cout << "ESC - exit\n";
}

//экспорт модели в VRML формат
void exportToVRML(const std::vector<double>& vx, const std::vector<double>& vy, const std::vector<double>& vz, const std::vector<int>& faces)
{
    writeVRML(vx, vy, vz, faces);
    std::cout << "Total vertices: " << vx.size() << std::endl;
    std::cout << "Total faces: " << faces.size() / 4 << std::endl;
}

//очистка
void cleanup(GLFWwindow* window)
{
    glfwTerminate();
}

int main()
{
    setlocale(LC_ALL, "Russian");

    //инициализация GLFW
    GLFWwindow* window = initGLFW();
    if (!window)
    {
        std::cout << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    printControls();

    //загрузка карты глубины
    std::string filePath = std::filesystem::current_path().string() + "\\resources\\DepthMap_7.dat";
    std::vector<std::vector<double>> dmap;
    int height, width;

    if (!loadDepthMap(filePath, dmap, height, width))
    {
        cleanup(window);
        return -1;
    }

    //преобразование в 3D координаты
    std::vector<double> vx, vy, vz;
    convertTo3DCoordinates(dmap, vx, vy, vz);
    std::vector<int> faces;
    createTriangleMesh(height, width, dmap, faces);

    //цикл рендеринга
    while (!glfwWindowShouldClose(window))
    {
        renderFrame(window, vx, vy, vz, faces);
    }

    //экспорт модели
    exportToVRML(vx, vy, vz, faces);

    cleanup(window);
    return 0;
}
