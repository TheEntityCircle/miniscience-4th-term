#include <iostream>
#include <cmath>
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

using namespace std;

// Класс расчётной точки
class CalcNode
{
// Класс сетки будет friend-ом точки
friend class CalcMesh;

protected:
    // Координаты
    double x;
    double y;
    double z;
    // Некая величина, в попугаях
    double smth;
    // Скорость
    double vx;
    double vy;
    double vz;

public:
    // Конструктор по умолчанию
    CalcNode() : x(0.0), y(0.0), z(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0)
    {
    }

    // Конструктор с указанием всех параметров
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz) 
            : x(x), y(y), z(z), smth(smth), vx(vx), vy(vy), vz(vz)
    {
    }

    // Метод отвечает за перемещение точки
    // Движемся время tau из текущего положения с текущей скоростью
    void move(double tau) {
        x += vx * tau;
        y += vy * tau;
        z += vz * tau;
    }
};

// Класс расчётной сетки
class CalcMesh
{
protected:
    // 3D-сетка из расчётных точек
    vector<vector<vector<CalcNode>>> points;

public:
    // Конструктор сетки size x size точек с шагом h по пространству
    CalcMesh(unsigned int size, double h) {
        points.resize(size);
        for(unsigned int i = 0; i < size; i++) {
            points[i].resize(size);
            for(unsigned int j = 0; j < size; j++) {
                points[i][j].resize(size);
                for(unsigned int k = 0; k < size; k++) {
                    // Начальные координаты зададим равномерно в XYZ
                    double pointX = i * h;
                    double pointY = j * h;
                    double pointZ = k * h;
                    // Модельная скалярная величина распределена как-то вот так
                    double smth = pow(pointX, 2) + pow(pointY, 2) + pow(pointZ, 2);
                    // Задавать скорости в данном примере поленимся, будут нули
                    points[i][j][k] = CalcNode(pointX, pointY, pointZ, smth, 0.0, 0.0, 0.0);
                }
            }
        }
    }

    // Метод отвечает за выполнение для всей сетки шага по времени величиной tau
    void doTimeStep(double tau) {
        // По сути метод просто двигает все точки
        for(unsigned int i = 0; i < points.size(); i++) {
            for(unsigned int j = 0; j < points[i].size(); j++) {
                for(unsigned int k = 0; k < points[i][j].size(); k++) {
                    points[i][j][k].move(tau);
                }
            }
        }
    }

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number) {
        // Сетка в терминах VTK
        vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        // Точки сетки в терминах VTK
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Скалярное поле на точках сетки
        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        // Векторное поле на точках сетки
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        // Обходим все точки нашей расчётной сетки
        unsigned int number = (unsigned int)points.size();
        for(unsigned int i = 0; i < number; i++) {
            for(unsigned int j = 0; j < number; j++) {
                for(unsigned int k = 0; k < number; k++) {
                    // Вставляем новую точку в сетку VTK-снапшота
                    dumpPoints->InsertNextPoint(points[i][j][k].x, points[i][j][k].y, points[i][j][k].z);

                    // Добавляем значение векторного поля в этой точке
                    double _vel[3] = {points[i][j][k].vx, points[i][j][k].vy, points[i][j][k].vz};
                    vel->InsertNextTuple(_vel);

                    // И значение скалярного поля тоже
                    smth->InsertNextValue(points[i][j][k].smth);
                }
            }
        }

        // Задаём размеры VTK-сетки (в точках, по трём осям)
        structuredGrid->SetDimensions(number, number, number);
        // Грузим точки в сетку
        structuredGrid->SetPoints(dumpPoints);

        // Присоединяем векторное и скалярное поля к точкам
        structuredGrid->GetPointData()->AddArray(vel);
        structuredGrid->GetPointData()->AddArray(smth);

        // Создаём снапшот в файле с заданным именем
        string fileName = "cubic3d-step-" + std::to_string(snap_number) + ".vts";
        vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(structuredGrid);
        writer->Write();
    }
};

int main()
{
    // Размер расчётной сетки, точек на сторону
    unsigned int size = 10;
    // Шаг точек по пространству
    double h = 0.1;
    // Шаг по времени
    double tau = 0.01;

    // Создаём сетку заданного размера
    CalcMesh mesh(size, h);

    // Пишем её начальное состояние в VTK
    mesh.snapshot(0);

    // Дальше можно сделать какие-нибудь шаги по времени аналогично 2D-примеру.
    // Но в этом примере их нет.

    return 0;
}
