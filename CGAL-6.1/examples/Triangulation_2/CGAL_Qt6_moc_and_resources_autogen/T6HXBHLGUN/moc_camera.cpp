/****************************************************************************
** Meta object code from reading C++ file 'camera.h'
**
** Created by: The Qt Meta Object Compiler version 68 (Qt 6.4.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../../../include/CGAL/Qt/camera.h"
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'camera.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 68
#error "This file was generated using the moc from 6.4.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

#ifndef Q_CONSTINIT
#define Q_CONSTINIT
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
namespace {
struct qt_meta_stringdata_CGAL__qglviewer__Camera_t {
    uint offsetsAndSizes[150];
    char stringdata0[24];
    char stringdata1[12];
    char stringdata2[1];
    char stringdata3[4];
    char stringdata4[4];
    char stringdata5[15];
    char stringdata6[11];
    char stringdata7[2];
    char stringdata8[6];
    char stringdata9[4];
    char stringdata10[12];
    char stringdata11[3];
    char stringdata12[7];
    char stringdata13[17];
    char stringdata14[10];
    char stringdata15[7];
    char stringdata16[7];
    char stringdata17[16];
    char stringdata18[10];
    char stringdata19[7];
    char stringdata20[7];
    char stringdata21[15];
    char stringdata22[4];
    char stringdata23[4];
    char stringdata24[16];
    char stringdata25[10];
    char stringdata26[12];
    char stringdata27[25];
    char stringdata28[6];
    char stringdata29[22];
    char stringdata30[14];
    char stringdata31[6];
    char stringdata32[3];
    char stringdata33[9];
    char stringdata34[8];
    char stringdata35[5];
    char stringdata36[5];
    char stringdata37[15];
    char stringdata38[4];
    char stringdata39[25];
    char stringdata40[5];
    char stringdata41[17];
    char stringdata42[15];
    char stringdata43[7];
    char stringdata44[24];
    char stringdata45[6];
    char stringdata46[7];
    char stringdata47[17];
    char stringdata48[20];
    char stringdata49[5];
    char stringdata50[24];
    char stringdata51[14];
    char stringdata52[2];
    char stringdata53[11];
    char stringdata54[15];
    char stringdata55[15];
    char stringdata56[24];
    char stringdata57[20];
    char stringdata58[14];
    char stringdata59[6];
    char stringdata60[23];
    char stringdata61[9];
    char stringdata62[24];
    char stringdata63[4];
    char stringdata64[24];
    char stringdata65[2];
    char stringdata66[22];
    char stringdata67[4];
    char stringdata68[18];
    char stringdata69[9];
    char stringdata70[11];
    char stringdata71[10];
    char stringdata72[12];
    char stringdata73[6];
    char stringdata74[16];
};
#define QT_MOC_LITERAL(ofs, len) \
    uint(sizeof(qt_meta_stringdata_CGAL__qglviewer__Camera_t::offsetsAndSizes) + ofs), len 
Q_CONSTINIT static const qt_meta_stringdata_CGAL__qglviewer__Camera_t qt_meta_stringdata_CGAL__qglviewer__Camera = {
    {
        QT_MOC_LITERAL(0, 23),  // "CGAL::qglviewer::Camera"
        QT_MOC_LITERAL(24, 11),  // "setPosition"
        QT_MOC_LITERAL(36, 0),  // ""
        QT_MOC_LITERAL(37, 3),  // "Vec"
        QT_MOC_LITERAL(41, 3),  // "pos"
        QT_MOC_LITERAL(45, 14),  // "setOrientation"
        QT_MOC_LITERAL(60, 10),  // "Quaternion"
        QT_MOC_LITERAL(71, 1),  // "q"
        QT_MOC_LITERAL(73, 5),  // "theta"
        QT_MOC_LITERAL(79, 3),  // "phi"
        QT_MOC_LITERAL(83, 11),  // "setUpVector"
        QT_MOC_LITERAL(95, 2),  // "up"
        QT_MOC_LITERAL(98, 6),  // "noMove"
        QT_MOC_LITERAL(105, 16),  // "setViewDirection"
        QT_MOC_LITERAL(122, 9),  // "direction"
        QT_MOC_LITERAL(132, 6),  // "lookAt"
        QT_MOC_LITERAL(139, 6),  // "target"
        QT_MOC_LITERAL(146, 15),  // "showEntireScene"
        QT_MOC_LITERAL(162, 9),  // "fitSphere"
        QT_MOC_LITERAL(172, 6),  // "center"
        QT_MOC_LITERAL(179, 6),  // "radius"
        QT_MOC_LITERAL(186, 14),  // "fitBoundingBox"
        QT_MOC_LITERAL(201, 3),  // "min"
        QT_MOC_LITERAL(205, 3),  // "max"
        QT_MOC_LITERAL(209, 15),  // "fitScreenRegion"
        QT_MOC_LITERAL(225, 9),  // "rectangle"
        QT_MOC_LITERAL(235, 11),  // "centerScene"
        QT_MOC_LITERAL(247, 24),  // "interpolateToZoomOnPixel"
        QT_MOC_LITERAL(272, 5),  // "pixel"
        QT_MOC_LITERAL(278, 21),  // "interpolateToFitScene"
        QT_MOC_LITERAL(300, 13),  // "interpolateTo"
        QT_MOC_LITERAL(314, 5),  // "Frame"
        QT_MOC_LITERAL(320, 2),  // "fr"
        QT_MOC_LITERAL(323, 8),  // "duration"
        QT_MOC_LITERAL(332, 7),  // "setType"
        QT_MOC_LITERAL(340, 4),  // "Type"
        QT_MOC_LITERAL(345, 4),  // "type"
        QT_MOC_LITERAL(350, 14),  // "setFieldOfView"
        QT_MOC_LITERAL(365, 3),  // "fov"
        QT_MOC_LITERAL(369, 24),  // "setHorizontalFieldOfView"
        QT_MOC_LITERAL(394, 4),  // "hfov"
        QT_MOC_LITERAL(399, 16),  // "setFOVToFitScene"
        QT_MOC_LITERAL(416, 14),  // "setAspectRatio"
        QT_MOC_LITERAL(431, 6),  // "aspect"
        QT_MOC_LITERAL(438, 23),  // "setScreenWidthAndHeight"
        QT_MOC_LITERAL(462, 5),  // "width"
        QT_MOC_LITERAL(468, 6),  // "height"
        QT_MOC_LITERAL(475, 16),  // "devicePixelRatio"
        QT_MOC_LITERAL(492, 19),  // "setZNearCoefficient"
        QT_MOC_LITERAL(512, 4),  // "coef"
        QT_MOC_LITERAL(517, 23),  // "setZClippingCoefficient"
        QT_MOC_LITERAL(541, 13),  // "setOrthoZNear"
        QT_MOC_LITERAL(555, 1),  // "z"
        QT_MOC_LITERAL(557, 10),  // "orthoZNear"
        QT_MOC_LITERAL(568, 14),  // "setSceneRadius"
        QT_MOC_LITERAL(583, 14),  // "setSceneCenter"
        QT_MOC_LITERAL(598, 23),  // "setSceneCenterFromPixel"
        QT_MOC_LITERAL(622, 19),  // "setSceneBoundingBox"
        QT_MOC_LITERAL(642, 13),  // "setPivotPoint"
        QT_MOC_LITERAL(656, 5),  // "point"
        QT_MOC_LITERAL(662, 22),  // "setPivotPointFromPixel"
        QT_MOC_LITERAL(685, 8),  // "setFrame"
        QT_MOC_LITERAL(694, 23),  // "ManipulatedCameraFrame*"
        QT_MOC_LITERAL(718, 3),  // "mcf"
        QT_MOC_LITERAL(722, 23),  // "setKeyFrameInterpolator"
        QT_MOC_LITERAL(746, 1),  // "i"
        QT_MOC_LITERAL(748, 21),  // "KeyFrameInterpolator*"
        QT_MOC_LITERAL(770, 3),  // "kfi"
        QT_MOC_LITERAL(774, 17),  // "addKeyFrameToPath"
        QT_MOC_LITERAL(792, 8),  // "playPath"
        QT_MOC_LITERAL(801, 10),  // "deletePath"
        QT_MOC_LITERAL(812, 9),  // "resetPath"
        QT_MOC_LITERAL(822, 11),  // "setFlySpeed"
        QT_MOC_LITERAL(834, 5),  // "speed"
        QT_MOC_LITERAL(840, 15)   // "onFrameModified"
    },
    "CGAL::qglviewer::Camera",
    "setPosition",
    "",
    "Vec",
    "pos",
    "setOrientation",
    "Quaternion",
    "q",
    "theta",
    "phi",
    "setUpVector",
    "up",
    "noMove",
    "setViewDirection",
    "direction",
    "lookAt",
    "target",
    "showEntireScene",
    "fitSphere",
    "center",
    "radius",
    "fitBoundingBox",
    "min",
    "max",
    "fitScreenRegion",
    "rectangle",
    "centerScene",
    "interpolateToZoomOnPixel",
    "pixel",
    "interpolateToFitScene",
    "interpolateTo",
    "Frame",
    "fr",
    "duration",
    "setType",
    "Type",
    "type",
    "setFieldOfView",
    "fov",
    "setHorizontalFieldOfView",
    "hfov",
    "setFOVToFitScene",
    "setAspectRatio",
    "aspect",
    "setScreenWidthAndHeight",
    "width",
    "height",
    "devicePixelRatio",
    "setZNearCoefficient",
    "coef",
    "setZClippingCoefficient",
    "setOrthoZNear",
    "z",
    "orthoZNear",
    "setSceneRadius",
    "setSceneCenter",
    "setSceneCenterFromPixel",
    "setSceneBoundingBox",
    "setPivotPoint",
    "point",
    "setPivotPointFromPixel",
    "setFrame",
    "ManipulatedCameraFrame*",
    "mcf",
    "setKeyFrameInterpolator",
    "i",
    "KeyFrameInterpolator*",
    "kfi",
    "addKeyFrameToPath",
    "playPath",
    "deletePath",
    "resetPath",
    "setFlySpeed",
    "speed",
    "onFrameModified"
};
#undef QT_MOC_LITERAL
} // unnamed namespace

Q_CONSTINIT static const uint qt_meta_data_CGAL__qglviewer__Camera[] = {

 // content:
      10,       // revision
       0,       // classname
       0,    0, // classinfo
      40,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags, initial metatype offsets
       1,    1,  254,    2, 0x0a,    1 /* Public */,
       5,    1,  257,    2, 0x0a,    3 /* Public */,
       5,    2,  260,    2, 0x0a,    5 /* Public */,
      10,    2,  265,    2, 0x0a,    8 /* Public */,
      10,    1,  270,    2, 0x2a,   11 /* Public | MethodCloned */,
      13,    1,  273,    2, 0x0a,   13 /* Public */,
      15,    1,  276,    2, 0x0a,   15 /* Public */,
      17,    0,  279,    2, 0x0a,   17 /* Public */,
      18,    2,  280,    2, 0x0a,   18 /* Public */,
      21,    2,  285,    2, 0x0a,   21 /* Public */,
      24,    1,  290,    2, 0x0a,   24 /* Public */,
      26,    0,  293,    2, 0x0a,   26 /* Public */,
      27,    1,  294,    2, 0x0a,   27 /* Public */,
      29,    0,  297,    2, 0x0a,   29 /* Public */,
      30,    2,  298,    2, 0x0a,   30 /* Public */,
      34,    1,  303,    2, 0x0a,   33 /* Public */,
      37,    1,  306,    2, 0x0a,   35 /* Public */,
      39,    1,  309,    2, 0x0a,   37 /* Public */,
      41,    0,  312,    2, 0x0a,   39 /* Public */,
      42,    1,  313,    2, 0x0a,   40 /* Public */,
      44,    3,  316,    2, 0x0a,   42 /* Public */,
      44,    2,  323,    2, 0x2a,   46 /* Public | MethodCloned */,
      48,    1,  328,    2, 0x0a,   49 /* Public */,
      50,    1,  331,    2, 0x0a,   51 /* Public */,
      51,    1,  334,    2, 0x0a,   53 /* Public */,
      53,    0,  337,    2, 0x0a,   55 /* Public */,
      54,    1,  338,    2, 0x0a,   56 /* Public */,
      55,    1,  341,    2, 0x0a,   58 /* Public */,
      56,    1,  344,    2, 0x0a,   60 /* Public */,
      57,    2,  347,    2, 0x0a,   62 /* Public */,
      58,    1,  352,    2, 0x0a,   65 /* Public */,
      60,    1,  355,    2, 0x0a,   67 /* Public */,
      61,    1,  358,    2, 0x0a,   69 /* Public */,
      64,    2,  361,    2, 0x0a,   71 /* Public */,
      68,    1,  366,    2, 0x0a,   74 /* Public */,
      69,    1,  369,    2, 0x0a,   76 /* Public */,
      70,    1,  372,    2, 0x0a,   78 /* Public */,
      71,    1,  375,    2, 0x0a,   80 /* Public */,
      72,    1,  378,    2, 0x0a,   82 /* Public */,
      74,    0,  381,    2, 0x08,   84 /* Private */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void, 0x80000000 | 6,    7,
    QMetaType::Void, QMetaType::QReal, QMetaType::QReal,    8,    9,
    QMetaType::Void, 0x80000000 | 3, QMetaType::Bool,   11,   12,
    QMetaType::Void, 0x80000000 | 3,   11,
    QMetaType::Void, 0x80000000 | 3,   14,
    QMetaType::Void, 0x80000000 | 3,   16,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 3, QMetaType::QReal,   19,   20,
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 3,   22,   23,
    QMetaType::Void, QMetaType::QRect,   25,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QPoint,   28,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 31, QMetaType::QReal,   32,   33,
    QMetaType::Void, 0x80000000 | 35,   36,
    QMetaType::Void, QMetaType::QReal,   38,
    QMetaType::Void, QMetaType::QReal,   40,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QReal,   43,
    QMetaType::Void, QMetaType::Int, QMetaType::Int, QMetaType::QReal,   45,   46,   47,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,   45,   46,
    QMetaType::Void, QMetaType::QReal,   49,
    QMetaType::Void, QMetaType::QReal,   49,
    QMetaType::Void, QMetaType::QReal,   52,
    QMetaType::QReal,
    QMetaType::Void, QMetaType::QReal,   20,
    QMetaType::Void, 0x80000000 | 3,   19,
    QMetaType::Bool, QMetaType::QPoint,   28,
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 3,   22,   23,
    QMetaType::Void, 0x80000000 | 3,   59,
    QMetaType::Bool, QMetaType::QPoint,   28,
    QMetaType::Void, 0x80000000 | 62,   63,
    QMetaType::Void, QMetaType::UInt, 0x80000000 | 66,   65,   67,
    QMetaType::Void, QMetaType::UInt,   65,
    QMetaType::Void, QMetaType::UInt,   65,
    QMetaType::Void, QMetaType::UInt,   65,
    QMetaType::Void, QMetaType::UInt,   65,
    QMetaType::Void, QMetaType::QReal,   73,
    QMetaType::Void,

       0        // eod
};

Q_CONSTINIT const QMetaObject CGAL::qglviewer::Camera::staticMetaObject = { {
    QMetaObject::SuperData::link<QObject::staticMetaObject>(),
    qt_meta_stringdata_CGAL__qglviewer__Camera.offsetsAndSizes,
    qt_meta_data_CGAL__qglviewer__Camera,
    qt_static_metacall,
    nullptr,
    qt_incomplete_metaTypeArray<qt_meta_stringdata_CGAL__qglviewer__Camera_t,
        // Q_OBJECT / Q_GADGET
        QtPrivate::TypeAndForceComplete<Camera, std::true_type>,
        // method 'setPosition'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        // method 'setOrientation'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Quaternion &, std::false_type>,
        // method 'setOrientation'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setUpVector'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        QtPrivate::TypeAndForceComplete<bool, std::false_type>,
        // method 'setUpVector'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        // method 'setViewDirection'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        // method 'lookAt'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        // method 'showEntireScene'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        // method 'fitSphere'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'fitBoundingBox'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        // method 'fitScreenRegion'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const QRect &, std::false_type>,
        // method 'centerScene'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        // method 'interpolateToZoomOnPixel'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const QPoint &, std::false_type>,
        // method 'interpolateToFitScene'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        // method 'interpolateTo'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Frame &, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setType'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<Type, std::false_type>,
        // method 'setFieldOfView'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setHorizontalFieldOfView'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setFOVToFitScene'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        // method 'setAspectRatio'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setScreenWidthAndHeight'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<int, std::false_type>,
        QtPrivate::TypeAndForceComplete<int, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setScreenWidthAndHeight'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<int, std::false_type>,
        QtPrivate::TypeAndForceComplete<int, std::false_type>,
        // method 'setZNearCoefficient'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setZClippingCoefficient'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setOrthoZNear'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'orthoZNear'
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setSceneRadius'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'setSceneCenter'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        // method 'setSceneCenterFromPixel'
        QtPrivate::TypeAndForceComplete<bool, std::false_type>,
        QtPrivate::TypeAndForceComplete<const QPoint &, std::false_type>,
        // method 'setSceneBoundingBox'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        // method 'setPivotPoint'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<const Vec &, std::false_type>,
        // method 'setPivotPointFromPixel'
        QtPrivate::TypeAndForceComplete<bool, std::false_type>,
        QtPrivate::TypeAndForceComplete<const QPoint &, std::false_type>,
        // method 'setFrame'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<ManipulatedCameraFrame * const, std::false_type>,
        // method 'setKeyFrameInterpolator'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<unsigned int, std::false_type>,
        QtPrivate::TypeAndForceComplete<KeyFrameInterpolator * const, std::false_type>,
        // method 'addKeyFrameToPath'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<unsigned int, std::false_type>,
        // method 'playPath'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<unsigned int, std::false_type>,
        // method 'deletePath'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<unsigned int, std::false_type>,
        // method 'resetPath'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<unsigned int, std::false_type>,
        // method 'setFlySpeed'
        QtPrivate::TypeAndForceComplete<void, std::false_type>,
        QtPrivate::TypeAndForceComplete<qreal, std::false_type>,
        // method 'onFrameModified'
        QtPrivate::TypeAndForceComplete<void, std::false_type>
    >,
    nullptr
} };

void CGAL::qglviewer::Camera::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Camera *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->setPosition((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1]))); break;
        case 1: _t->setOrientation((*reinterpret_cast< std::add_pointer_t<Quaternion>>(_a[1]))); break;
        case 2: _t->setOrientation((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<qreal>>(_a[2]))); break;
        case 3: _t->setUpVector((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<bool>>(_a[2]))); break;
        case 4: _t->setUpVector((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1]))); break;
        case 5: _t->setViewDirection((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1]))); break;
        case 6: _t->lookAt((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1]))); break;
        case 7: _t->showEntireScene(); break;
        case 8: _t->fitSphere((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<qreal>>(_a[2]))); break;
        case 9: _t->fitBoundingBox((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<Vec>>(_a[2]))); break;
        case 10: _t->fitScreenRegion((*reinterpret_cast< std::add_pointer_t<QRect>>(_a[1]))); break;
        case 11: _t->centerScene(); break;
        case 12: _t->interpolateToZoomOnPixel((*reinterpret_cast< std::add_pointer_t<QPoint>>(_a[1]))); break;
        case 13: _t->interpolateToFitScene(); break;
        case 14: _t->interpolateTo((*reinterpret_cast< std::add_pointer_t<Frame>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<qreal>>(_a[2]))); break;
        case 15: _t->setType((*reinterpret_cast< std::add_pointer_t<Type>>(_a[1]))); break;
        case 16: _t->setFieldOfView((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1]))); break;
        case 17: _t->setHorizontalFieldOfView((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1]))); break;
        case 18: _t->setFOVToFitScene(); break;
        case 19: _t->setAspectRatio((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1]))); break;
        case 20: _t->setScreenWidthAndHeight((*reinterpret_cast< std::add_pointer_t<int>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<int>>(_a[2])),(*reinterpret_cast< std::add_pointer_t<qreal>>(_a[3]))); break;
        case 21: _t->setScreenWidthAndHeight((*reinterpret_cast< std::add_pointer_t<int>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<int>>(_a[2]))); break;
        case 22: _t->setZNearCoefficient((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1]))); break;
        case 23: _t->setZClippingCoefficient((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1]))); break;
        case 24: _t->setOrthoZNear((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1]))); break;
        case 25: { qreal _r = _t->orthoZNear();
            if (_a[0]) *reinterpret_cast< qreal*>(_a[0]) = std::move(_r); }  break;
        case 26: _t->setSceneRadius((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1]))); break;
        case 27: _t->setSceneCenter((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1]))); break;
        case 28: { bool _r = _t->setSceneCenterFromPixel((*reinterpret_cast< std::add_pointer_t<QPoint>>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = std::move(_r); }  break;
        case 29: _t->setSceneBoundingBox((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<Vec>>(_a[2]))); break;
        case 30: _t->setPivotPoint((*reinterpret_cast< std::add_pointer_t<Vec>>(_a[1]))); break;
        case 31: { bool _r = _t->setPivotPointFromPixel((*reinterpret_cast< std::add_pointer_t<QPoint>>(_a[1])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = std::move(_r); }  break;
        case 32: _t->setFrame((*reinterpret_cast< std::add_pointer_t<ManipulatedCameraFrame*>>(_a[1]))); break;
        case 33: _t->setKeyFrameInterpolator((*reinterpret_cast< std::add_pointer_t<uint>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<KeyFrameInterpolator*>>(_a[2]))); break;
        case 34: _t->addKeyFrameToPath((*reinterpret_cast< std::add_pointer_t<uint>>(_a[1]))); break;
        case 35: _t->playPath((*reinterpret_cast< std::add_pointer_t<uint>>(_a[1]))); break;
        case 36: _t->deletePath((*reinterpret_cast< std::add_pointer_t<uint>>(_a[1]))); break;
        case 37: _t->resetPath((*reinterpret_cast< std::add_pointer_t<uint>>(_a[1]))); break;
        case 38: _t->setFlySpeed((*reinterpret_cast< std::add_pointer_t<qreal>>(_a[1]))); break;
        case 39: _t->onFrameModified(); break;
        default: ;
        }
    }
}

const QMetaObject *CGAL::qglviewer::Camera::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *CGAL::qglviewer::Camera::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_CGAL__qglviewer__Camera.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int CGAL::qglviewer::Camera::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 40)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 40;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 40)
            *reinterpret_cast<QMetaType *>(_a[0]) = QMetaType();
        _id -= 40;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
