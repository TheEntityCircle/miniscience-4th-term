
#include <set>
#include <iostream>
#include <gmsh.h>
#include <vector>



int main(int argc, char **argv)
{

    gmsh::initialize();
    gmsh::model::add("cilinder_model");
    double MAINCenterX = 0.0; // X- координата центра большого круга
    double MAINCenterY = 0.0; // Y- координата центра большого круга
    double MAINCenterZ = 0.0; // Z - координата центра большого круга
    double MAINradius = 1.0; // Радиус большого круга
    
    double LEFTCenterX = MAINCenterX ; // X- координата центра левого круга
    double LEFTCenterY = MAINCenterY - MAINradius; // Y- координата центра левого круга
    double LEFTCenterZ = MAINCenterZ; // Z - координата центра левого круга
    
    
    double RIGHTCenterX = MAINCenterX ; // X- координата центра правого круга
    double RIGHTCenterY = MAINCenterY + MAINradius; // Y- координата центра правого круга
    double RIGHTCenterZ = MAINCenterZ; // Z - координата центра правого круга
    
    
    double radius = 0.1; // Радиус меньшего круга
    
    
    
    double lc = 0.01;
    
    // Создание точек для верхнего круга
    int centerTag_UP = gmsh::model::geo::addPoint(MAINCenterX, MAINCenterY, MAINCenterZ + radius, lc); //1 Точка
    int startTag_UP = gmsh::model::geo::addPoint(MAINCenterX, MAINCenterY + MAINradius, MAINCenterZ + radius, lc); //2 Точка справа сверху
    int endTag_UP = gmsh::model::geo::addPoint(MAINCenterX , MAINCenterY - MAINradius , MAINCenterZ + radius, lc); //3 Точка слева сверху
    
    //int circleTag_UP1 = gmsh::model::geo::addCircleArc(startTag_UP, centerTag_UP, endTag_UP); //1
    int circleTag_UP2 = gmsh::model::geo::addCircleArc(endTag_UP, centerTag_UP, startTag_UP); //2
    
    //int curveLoopTag_UP = gmsh::model::geo::addCurveLoop({circleTag_UP1, circleTag_UP2});
    
    //Создаем нижнего круга
    
    int centerTag_DOWN = gmsh::model::geo::addPoint(MAINCenterX, MAINCenterY, MAINCenterZ - radius, lc); // 4 точка
    int startTag_DOWN = gmsh::model::geo::addPoint(MAINCenterX , MAINCenterY + MAINradius, MAINCenterZ - radius, lc); // 5 Справа снизу
    int endTag_DOWN = gmsh::model::geo::addPoint(MAINCenterX , MAINCenterY  - MAINradius, MAINCenterZ - radius, lc); // 6 Слева снизу
    
    //int circleTag_DOWN1 = gmsh::model::geo::addCircleArc(startTag_DOWN, centerTag_DOWN, endTag_DOWN); //3
    int circleTag_DOWN2 = gmsh::model::geo::addCircleArc(endTag_DOWN, centerTag_DOWN, startTag_DOWN); //4
    
    //int curveLoopunderTag_DOWN = gmsh::model::geo::addCurveLoop({circleTag_DOWN1, circleTag_DOWN2});
   
    
    
    
    
    
    //Создаем левый круг
    
    int centerTag_left = gmsh::model::geo::addPoint(LEFTCenterX, LEFTCenterY, LEFTCenterZ, lc); // 7 точка
    int upTag_left = endTag_UP; // 3 точка(была 8)
    int innTag_left = gmsh::model::geo::addPoint(LEFTCenterX , LEFTCenterY + radius, LEFTCenterZ, lc); // 8 точка
    int downTag_left = endTag_DOWN; // 6 точка
    int extTag_left = gmsh::model::geo::addPoint(LEFTCenterX , LEFTCenterY - radius, LEFTCenterZ, lc);  // 9 точка
    
    
    
   // int circleTag_left1 = gmsh::model::geo::addCircleArc(upTag_left, centerTag_left, innTag_left); //5
   // int circleTag_left2 = gmsh::model::geo::addCircleArc(innTag_left, centerTag_left, downTag_left); //6
    int circleTag_left3 = gmsh::model::geo::addCircleArc(downTag_left, centerTag_left, extTag_left); //7
    int circleTag_left4 = gmsh::model::geo::addCircleArc(extTag_left, centerTag_left, upTag_left);
    
    
    
    //int CircleTag_inn_LEFT = gmsh::model::geo::addCompoundSpline({circleTag_left1, circleTag_left2});
   // int CircleTag_ext_LEFT = gmsh::model::geo::addCompoundSpline({circleTag_left3, circleTag_left4});
    
    
   // int curveLoopunderTag_left = gmsh::model::geo::addCurveLoop({circleTag_left1, circleTag_left2,circleTag_left3, circleTag_left4 });
    
    // Создаем правый круг
    
        int centerTag_RIGHT = gmsh::model::geo::addPoint(RIGHTCenterX, RIGHTCenterY, RIGHTCenterZ, lc); // 10 точка
        int upTag_RIGHT = startTag_UP; //2 точка
        int innTag_RIGHT = gmsh::model::geo::addPoint(RIGHTCenterX , RIGHTCenterY + radius, RIGHTCenterZ, lc); // 11 точка
    int downTag_RIGHT = startTag_DOWN; // 5 точка
        int extTag_RIGHT = gmsh::model::geo::addPoint(RIGHTCenterX , RIGHTCenterY - radius, RIGHTCenterZ, lc); // 12 точка
        
        
        
        int circleTag_RIGHT1 = gmsh::model::geo::addCircleArc(upTag_RIGHT, centerTag_RIGHT, innTag_RIGHT); //9
        int circleTag_RIGHT2 = gmsh::model::geo::addCircleArc(innTag_RIGHT, centerTag_RIGHT, downTag_RIGHT); //10
        //int circleTag_RIGHT3 = gmsh::model::geo::addCircleArc(downTag_RIGHT, centerTag_RIGHT, extTag_RIGHT); //11
        //int circleTag_RIGHT4 = gmsh::model::geo::addCircleArc(extTag_RIGHT, centerTag_RIGHT, upTag_RIGHT); // 12
    
    
    // Создаем внешний круг
    
    int centerTag_EXT = gmsh::model::geo::addPoint(MAINCenterX, MAINCenterY, MAINCenterZ, lc);
    int rightTag_EXT = innTag_RIGHT;
    int leftTag_EXT = extTag_left;
    int circleTag_EXT1 = gmsh::model::geo::addCircleArc(leftTag_EXT, centerTag_EXT, rightTag_EXT);
    
    
    //int CircleTag_ext_RIGHT = gmsh::model::geo::addCompoundSpline({circleTag_RIGHT1, circleTag_RIGHT2});
    //int CircleTag_inn_RIGHT = gmsh::model::geo::addCompoundSpline({circleTag_RIGHT3, circleTag_RIGHT4});
        
        
    //int curveLoopTag_RIGHT = gmsh::model::geo::addCurveLoop({circleTag_RIGHT1, circleTag_RIGHT2,circleTag_RIGHT3, circleTag_RIGHT4 });
    
    //int curveLoopTag_EXT = gmsh::model::geo::addCurveLoop({circleTag_UP2,circleTag_RIGHT1, circleTag_RIGHT2 ,-circleTag_DOWN2,- circleTag_left3, -circleTag_left4 });
    
    
    
    
    std::cout << "circleTag_UP2::" << circleTag_UP2 << "\n";
    std::cout << "circleTag_RIGHT1::" << circleTag_RIGHT1 << "\n";
    std::cout << "circleTag_RIGHT2::" << circleTag_RIGHT2 << "\n";
    std::cout << "circleTag_DOWN2::" << circleTag_DOWN2 << "\n";
    std::cout << "circleTag_left3::" << circleTag_left3 << "\n";
    std::cout << "circleTag_left4::" << circleTag_left4 << "\n";

    std::cout << "======================" << "\n";
    //std::vector<std::pair<int, int>> adj;
    //gmsh::model::getAdjacencies(1, circleTag_UP2, adj, false);
    //std::cout << "StartcircleTag_UP2::" << adj[0].second << "\n";
    //std::cout << "EndcircleTag_UP2::" << adj[1].second << "\n";
    
    int curveLoopTag_EXT_UP = gmsh::model::geo::addCurveLoop({circleTag_UP2, circleTag_RIGHT1,-circleTag_EXT1 , circleTag_left4});
    
    int surface = gmsh::model::geo::addPlaneSurface({curveLoopTag_EXT_UP});
    
    //int curveLoopTag_EXT = gmsh::model::geo::addCurveLoop({circleTag_UP2, circleTag_RIGHT1, circleTag_RIGHT2, -circleTag_DOWN2, circleTag_left3, circleTag_left4});
    //int surface = gmsh::model::geo::addPlaneSurface({curveLoopTag_EXT});
    //gmsh::model::geo::rotate({{1, circleTag_left1}, {1, circleTag_left2}, {0, centerTag_left}, {0, startTag_left}, {0, endTag_left}}, LEFTCenterX, LEFTCenterY, LEFTCenterZ, 0, 1, 0, angle);
    
    
    //Создаем правый круг
    /*
    int centerTag_right = gmsh::model::geo::addPoint(RIGHTCenterX, RIGHTCenterY, RIGHTCenterZ, lc);
    int startTag_right = gmsh::model::geo::addPoint(RIGHTCenterX , RIGHTCenterY - radius, RIGHTCenterZ, lc); //ближняя к центру(?)
    int endTag_right = gmsh::model::geo::addPoint(RIGHTCenterX , RIGHTCenterY + radius, MAINCenterZ, lc); // Дальная от центра(?)
    
    int circleTag_right1 = gmsh::model::geo::addCircleArc(startTag_right, centerTag_right, endTag_right); //7
    int circleTag_right2 = gmsh::model::geo::addCircleArc(endTag_right, centerTag_right, startTag_right); //8
    
    
    
    int curveLoopunderTag_right = gmsh::model::geo::addCurveLoop({circleTag_right1, circleTag_right2});
    
    */
    
    //Создаем верхний круг
    /*
    int centerupTag = gmsh::model::geo::addPoint(centerX, centerY, centerZup, lc);
    int startupTag = gmsh::model::geo::addPoint(centerX + radius, centerY, centerZup, lc);
    int endupTag = gmsh::model::geo::addPoint(centerX - radius, centerY, centerZup, lc);
    
    int circleupTag1 = gmsh::model::geo::addCircleArc(startupTag, centerupTag, endupTag); //3
    int circleupTag2 = gmsh::model::geo::addCircleArc(endupTag, centerupTag, startupTag); //4
    
    int curveLoopupTag = gmsh::model::geo::addCurveLoop({circleupTag1, circleupTag2});
    
    // Создаем линии для поверхности посередине
    int firstmedlineTag = gmsh::model::geo::addLine(startunderTag, startupTag); //5
    
    int secondmedlineTag = gmsh::model::geo::addLine(endunderTag, endupTag); //6
    
    int curveLoopMedTag1 = gmsh::model::geo::addCurveLoop({circleupTag1, -secondmedlineTag, -circleunderTag1, firstmedlineTag});
    
    int curveLoopMedTag2 = gmsh::model::geo::addCurveLoop({circleupTag2, -firstmedlineTag, -circleunderTag2, secondmedlineTag});

    // Создание поверхностей
    
    
    // Создаем нижний круг
    int surfaceunderTag = gmsh::model::geo::addPlaneSurface({curveLoopunderTag});
    //Создаем верхний круг
    int surfaceupTag = gmsh::model::geo::addPlaneSurface({curveLoopupTag});
    
    // Создаем половину середины
    int surfacemediumTag1 = gmsh::model::geo::addSurfaceFilling({curveLoopMedTag1});
    // И вторую половину
    int surfacemediumTag2 = gmsh::model::geo::addSurfaceFilling({curveLoopMedTag2});
    
    //Остался Объем
    int SurfaceLoopTag1 = gmsh::model::geo::addSurfaceLoop({surfaceunderTag, surfacemediumTag1, surfacemediumTag2, surfaceupTag});
    
    int ValueTag = gmsh::model::geo::addVolume({SurfaceLoopTag1});
     */
    // Синхронизация геометрии
    gmsh::model::geo::synchronize();

    // Генерация 2D сетки
    gmsh::model::mesh::generate(2);

  
    gmsh::write("circle_model.msh");

    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) gmsh::fltk::run();

  
    gmsh::finalize();

    return 0;
}




