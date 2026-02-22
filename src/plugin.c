/*  plugin.c  —  jp-physics v1.0.0
 *  はじむ言語プラグインバインディング (Unity レベル物理エンジン)
 */
#include "../../jp/include/hajimu_plugin.h"
#include "../include/phys.h"
#include <string.h>
#include <stdlib.h>

/* ─── ヘルパーマクロ ─────────────────────────────────── */
static PhysWorld* g_world = NULL;

static double NUM(Value* v){ return (v&&v->type==VALUE_NUMBER)?v->number:0.0; }
static int    INT(Value* v){ return (int)NUM(v); }
static int    BOL_V(Value* v){ return (v&&v->type==VALUE_BOOL)?(int)v->boolean:(v&&v->type==VALUE_NUMBER)?v->number!=0:0; }
#define vN(n)  hajimu_number(n)
#define vB(b)  hajimu_bool((bool)(b))
#define vNULL() hajimu_null()

static Value vArray3(float x, float y, float z){
    Value arr = hajimu_array();
    hajimu_array_push(&arr, hajimu_number(x));
    hajimu_array_push(&arr, hajimu_number(y));
    hajimu_array_push(&arr, hajimu_number(z));
    return arr;
}

/* ═══════════════════════════════════════════════════════
 *  ワールド
 * ═══════════════════════════════════════════════════════ */
static Value p_world_create(int argc, Value* argv){
    (void)argc;(void)argv;
    if (g_world) phys_world_destroy(g_world);
    g_world = phys_world_create();
    return vB(g_world != NULL);
}
static Value p_world_destroy(int argc, Value* argv){
    (void)argc;(void)argv;
    if (g_world){ phys_world_destroy(g_world); g_world=NULL; }
    return vNULL();
}
static Value p_world_step(int argc, Value* argv){
    float dt = argc>=1?(float)NUM(&argv[0]):(1.0f/60.0f);
    if (g_world) phys_world_step(g_world, dt);
    return vNULL();
}
static Value p_set_gravity(int argc, Value* argv){
    float x=argc>=1?(float)NUM(&argv[0]):0;
    float y=argc>=2?(float)NUM(&argv[1]):-9.81f;
    float z=argc>=3?(float)NUM(&argv[2]):0;
    if (g_world) phys_set_gravity(g_world, x, y, z);
    return vNULL();
}
static Value p_set_substeps(int argc, Value* argv){
    if (g_world&&argc>=1) phys_set_substeps(g_world, INT(&argv[0]));
    return vNULL();
}
static Value p_set_layer_col(int argc, Value* argv){
    if (!g_world||argc<3) return vNULL();
    phys_set_layer_collision(g_world, INT(&argv[0]), INT(&argv[1]), BOL_V(&argv[2]));
    return vNULL();
}

/* ═══════════════════════════════════════════════════════
 *  剛体生成 / 破棄
 * ═══════════════════════════════════════════════════════ */
static Value p_body_create(int argc, Value* argv){
    (void)argc;(void)argv;
    return vN(phys_body_create(g_world));
}
static Value p_body_destroy(int argc, Value* argv){
    if (g_world&&argc>=1) phys_body_destroy(g_world, INT(&argv[0]));
    return vNULL();
}

/* ─── 位置 / 回転 ─────────────────────────────────────── */
static Value p_set_pos(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_body_set_pos(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}
static Value p_get_pos(int argc, Value* argv){
    float x=0,y=0,z=0;
    if (g_world&&argc>=1) phys_body_get_pos(g_world,INT(&argv[0]),&x,&y,&z);
    return vArray3(x,y,z);
}
static Value p_set_rot(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_body_set_rot_euler(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}
static Value p_get_rot(int argc, Value* argv){
    float rx=0,ry=0,rz=0;
    if (g_world&&argc>=1) phys_body_get_rot_euler(g_world,INT(&argv[0]),&rx,&ry,&rz);
    return vArray3(rx,ry,rz);
}

/* ─── 速度 ────────────────────────────────────────────── */
static Value p_set_vel(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_body_set_vel(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}
static Value p_get_vel(int argc, Value* argv){
    float x=0,y=0,z=0;
    if (g_world&&argc>=1) phys_body_get_vel(g_world,INT(&argv[0]),&x,&y,&z);
    return vArray3(x,y,z);
}
static Value p_set_omega(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_body_set_omega(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}
static Value p_get_omega(int argc, Value* argv){
    float x=0,y=0,z=0;
    if (g_world&&argc>=1) phys_body_get_omega(g_world,INT(&argv[0]),&x,&y,&z);
    return vArray3(x,y,z);
}

/* ─── 質量 / マテリアル / 減衰 ───────────────────────── */
static Value p_set_mass(int argc, Value* argv){
    if (!g_world||argc<2) return vNULL();
    phys_body_set_mass(g_world,INT(&argv[0]),(float)NUM(&argv[1]));
    return vNULL();
}
static Value p_set_drag(int argc, Value* argv){
    if (!g_world||argc<3) return vNULL();
    phys_body_set_drag(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]));
    return vNULL();
}
static Value p_set_material(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_body_set_material(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}

/* ─── タイプ / フラグ ─────────────────────────────────── */
static Value p_set_type(int argc, Value* argv){
    if (!g_world||argc<2) return vNULL();
    phys_body_set_type(g_world,INT(&argv[0]),(PhysBodyType)INT(&argv[1]));
    return vNULL();
}
static Value p_set_trigger(int argc, Value* argv){
    if (!g_world||argc<2) return vNULL();
    phys_body_set_trigger(g_world,INT(&argv[0]),BOL_V(&argv[1]));
    return vNULL();
}
static Value p_set_layer(int argc, Value* argv){
    if (!g_world||argc<3) return vNULL();
    phys_body_set_layer(g_world,INT(&argv[0]),INT(&argv[1]),(uint32_t)INT(&argv[2]));
    return vNULL();
}
static Value p_freeze_pos(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_body_freeze_pos(g_world,INT(&argv[0]),BOL_V(&argv[1]),BOL_V(&argv[2]),BOL_V(&argv[3]));
    return vNULL();
}
static Value p_freeze_rot(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_body_freeze_rot(g_world,INT(&argv[0]),BOL_V(&argv[1]),BOL_V(&argv[2]),BOL_V(&argv[3]));
    return vNULL();
}
static Value p_set_sleep(int argc, Value* argv){
    if (!g_world||argc<2) return vNULL();
    phys_body_set_sleep(g_world,INT(&argv[0]),BOL_V(&argv[1]));
    return vNULL();
}
static Value p_is_sleeping(int argc, Value* argv){
    if (!g_world||argc<1) return vB(0);
    return vB(phys_body_is_sleeping(g_world,INT(&argv[0])));
}
static Value p_set_tag(int argc, Value* argv){
    if (!g_world||argc<2) return vNULL();
    phys_body_set_tag(g_world,INT(&argv[0]),INT(&argv[1]));
    return vNULL();
}
static Value p_get_tag(int argc, Value* argv){
    if (!g_world||argc<1) return vN(-1);
    return vN(phys_body_get_tag(g_world,INT(&argv[0])));
}

/* ─── コライダー ──────────────────────────────────────── */
static Value p_col_sphere(int argc, Value* argv){
    if (!g_world||argc<2) return vNULL();
    phys_col_sphere(g_world,INT(&argv[0]),(float)NUM(&argv[1]));
    return vNULL();
}
static Value p_col_box(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_col_box(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}
static Value p_col_capsule(int argc, Value* argv){
    if (!g_world||argc<3) return vNULL();
    phys_col_capsule(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]));
    return vNULL();
}
static Value p_col_offset(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_col_offset(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}

/* ─── 力 ──────────────────────────────────────────────── */
static Value p_add_force(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_add_force(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}
static Value p_add_impulse(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_add_impulse(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}
static Value p_add_torque(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    phys_add_torque(g_world,INT(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]));
    return vNULL();
}
static Value p_add_force_at(int argc, Value* argv){
    if (!g_world||argc<7) return vNULL();
    phys_add_force_at(g_world,INT(&argv[0]),
        (float)NUM(&argv[1]),(float)NUM(&argv[2]),(float)NUM(&argv[3]),
        (float)NUM(&argv[4]),(float)NUM(&argv[5]),(float)NUM(&argv[6]));
    return vNULL();
}
static Value p_clear_forces(int argc, Value* argv){
    if (!g_world||argc<1) return vNULL();
    phys_clear_forces(g_world,INT(&argv[0]));
    return vNULL();
}

/* ─── ジョイント ──────────────────────────────────────── */
static Value p_joint_fixed(int argc, Value* argv){
    if (!g_world||argc<2) return vN(-1);
    return vN(phys_joint_fixed(g_world,INT(&argv[0]),INT(&argv[1])));
}
static Value p_joint_distance(int argc, Value* argv){
    if (!g_world||argc<4) return vN(-1);
    return vN(phys_joint_distance(g_world,INT(&argv[0]),INT(&argv[1]),
                                   (float)NUM(&argv[2]),(float)NUM(&argv[3])));
}
static Value p_joint_spring(int argc, Value* argv){
    if (!g_world||argc<5) return vN(-1);
    return vN(phys_joint_spring(g_world,INT(&argv[0]),INT(&argv[1]),
                                  (float)NUM(&argv[2]),(float)NUM(&argv[3]),(float)NUM(&argv[4])));
}
static Value p_joint_hinge(int argc, Value* argv){
    if (!g_world||argc<7) return vN(-1);
    return vN(phys_joint_hinge(g_world,INT(&argv[0]),INT(&argv[1]),
                                 (float)NUM(&argv[2]),(float)NUM(&argv[3]),(float)NUM(&argv[4]),
                                 (float)NUM(&argv[5]),(float)NUM(&argv[6])));
}
static Value p_joint_destroy(int argc, Value* argv){
    if (!g_world||argc<1) return vNULL();
    phys_joint_destroy(g_world,INT(&argv[0]));
    return vNULL();
}

/* ─── レイキャスト / クエリ ───────────────────────────── */
static Value p_raycast(int argc, Value* argv){
    if (!g_world||argc<6) return vNULL();
    float max_d = argc>=7?(float)NUM(&argv[6]):1000.0f;
    uint32_t mask = argc>=8?(uint32_t)INT(&argv[7]):0xFFFFFFFF;
    PhysRayHit h = phys_raycast(g_world,
        (float)NUM(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),
        (float)NUM(&argv[3]),(float)NUM(&argv[4]),(float)NUM(&argv[5]),
        max_d, mask);
    /* 戻り値: [hit, body_id, t, nx, ny, nz] */
    Value arr = hajimu_array();
    hajimu_array_push(&arr, hajimu_bool(h.hit));
    hajimu_array_push(&arr, hajimu_number(h.body_id));
    hajimu_array_push(&arr, hajimu_number(h.t));
    hajimu_array_push(&arr, hajimu_number(h.normal.x));
    hajimu_array_push(&arr, hajimu_number(h.normal.y));
    hajimu_array_push(&arr, hajimu_number(h.normal.z));
    return arr;
}
static Value p_sphere_cast(int argc, Value* argv){
    if (!g_world||argc<7) return vNULL();
    float max_d = argc>=8?(float)NUM(&argv[7]):1000.0f;
    uint32_t mask = argc>=9?(uint32_t)INT(&argv[8]):0xFFFFFFFF;
    PhysRayHit h = phys_sphere_cast(g_world,
        (float)NUM(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),
        (float)NUM(&argv[3]),(float)NUM(&argv[4]),(float)NUM(&argv[5]),
        (float)NUM(&argv[6]), max_d, mask);
    Value arr = hajimu_array();
    hajimu_array_push(&arr, hajimu_bool(h.hit));
    hajimu_array_push(&arr, hajimu_number(h.body_id));
    hajimu_array_push(&arr, hajimu_number(h.t));
    hajimu_array_push(&arr, hajimu_number(h.normal.x));
    hajimu_array_push(&arr, hajimu_number(h.normal.y));
    hajimu_array_push(&arr, hajimu_number(h.normal.z));
    return arr;
}
static Value p_overlap_sphere(int argc, Value* argv){
    if (!g_world||argc<4) return vNULL();
    uint32_t mask = argc>=5?(uint32_t)INT(&argv[4]):0xFFFFFFFF;
    int results[PHYS_MAX_OVERLAP];
    int n = phys_overlap_sphere(g_world,
        (float)NUM(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),
        (float)NUM(&argv[3]), mask, results, PHYS_MAX_OVERLAP);
    Value arr = hajimu_array();
    for (int i=0;i<n;i++) hajimu_array_push(&arr, hajimu_number(results[i]));
    return arr;
}
static Value p_overlap_box(int argc, Value* argv){
    if (!g_world||argc<6) return vNULL();
    uint32_t mask = argc>=7?(uint32_t)INT(&argv[6]):0xFFFFFFFF;
    int results[PHYS_MAX_OVERLAP];
    int n = phys_overlap_box(g_world,
        (float)NUM(&argv[0]),(float)NUM(&argv[1]),(float)NUM(&argv[2]),
        (float)NUM(&argv[3]),(float)NUM(&argv[4]),(float)NUM(&argv[5]),
        mask, results, PHYS_MAX_OVERLAP);
    Value arr = hajimu_array();
    for (int i=0;i<n;i++) hajimu_array_push(&arr, hajimu_number(results[i]));
    return arr;
}
static Value p_check_col(int argc, Value* argv){
    if (!g_world||argc<2) return vB(0);
    return vB(phys_check_collision(g_world,INT(&argv[0]),INT(&argv[1])));
}
static Value p_check_trig(int argc, Value* argv){
    if (!g_world||argc<2) return vB(0);
    return vB(phys_check_trigger(g_world,INT(&argv[0]),INT(&argv[1])));
}

/* ─── ユーティリティ ──────────────────────────────────── */
static Value p_get_forward(int argc, Value* argv){
    float x=0,y=0,z=-1;
    if (g_world&&argc>=1) phys_body_get_forward(g_world,INT(&argv[0]),&x,&y,&z);
    return vArray3(x,y,z);
}
static Value p_get_up(int argc, Value* argv){
    float x=0,y=1,z=0;
    if (g_world&&argc>=1) phys_body_get_up(g_world,INT(&argv[0]),&x,&y,&z);
    return vArray3(x,y,z);
}
static Value p_get_right(int argc, Value* argv){
    float x=1,y=0,z=0;
    if (g_world&&argc>=1) phys_body_get_right(g_world,INT(&argv[0]),&x,&y,&z);
    return vArray3(x,y,z);
}
static Value p_get_speed(int argc, Value* argv){
    if (!g_world||argc<1) return vN(0);
    return vN(phys_body_get_speed(g_world,INT(&argv[0])));
}
static Value p_contact_count(int argc, Value* argv){
    (void)argc;(void)argv;
    return vN(phys_get_contact_count(g_world));
}
static Value p_contact_bodies(int argc, Value* argv){
    if (!g_world||argc<1) return vNULL();
    int a=-1,b=-1;
    phys_get_contact_bodies(g_world,INT(&argv[0]),&a,&b);
    Value arr = hajimu_array();
    hajimu_array_push(&arr, hajimu_number(a));
    hajimu_array_push(&arr, hajimu_number(b));
    return arr;
}
static Value p_contact_info(int argc, Value* argv){
    if (!g_world||argc<1) return vNULL();
    float nx=0,ny=1,nz=0,depth=0;
    phys_get_contact_info(g_world,INT(&argv[0]),&nx,&ny,&nz,&depth);
    Value arr = hajimu_array();
    hajimu_array_push(&arr, hajimu_number(nx));
    hajimu_array_push(&arr, hajimu_number(ny));
    hajimu_array_push(&arr, hajimu_number(nz));
    hajimu_array_push(&arr, hajimu_number(depth));
    return arr;
}

/* ═══════════════════════════════════════════════════════
 *  関数テーブル
 * ═══════════════════════════════════════════════════════ */
static HajimuPluginFunc functions[] = {
    /* ワールド */
    {"物理ワールド作成",   p_world_create,   0, 0},
    {"物理ワールド削除",   p_world_destroy,  0, 0},
    {"物理更新",          p_world_step,     0, 1},
    {"重力設定",          p_set_gravity,    3, 3},
    {"サブステップ設定",   p_set_substeps,   1, 1},
    {"レイヤー衝突設定",   p_set_layer_col,  3, 3},
    /* 剛体 */
    {"剛体作成",          p_body_create,    0, 0},
    {"剛体削除",          p_body_destroy,   1, 1},
    /* 位置 / 回転 */
    {"位置設定",          p_set_pos,        4, 4},
    {"位置取得",          p_get_pos,        1, 1},
    {"回転設定",          p_set_rot,        4, 4},
    {"回転取得",          p_get_rot,        1, 1},
    /* 速度 / 角速度 */
    {"速度設定",          p_set_vel,        4, 4},
    {"速度取得",          p_get_vel,        1, 1},
    {"角速度設定",        p_set_omega,      4, 4},
    {"角速度取得",        p_get_omega,      1, 1},
    /* 質量 / マテリアル */
    {"質量設定",          p_set_mass,       2, 2},
    {"減衰設定",          p_set_drag,       3, 3},
    {"物理マテリアル設定", p_set_material,   4, 4},
    /* タイプ / フラグ */
    {"剛体タイプ設定",    p_set_type,       2, 2},
    {"トリガー設定",      p_set_trigger,    2, 2},
    {"レイヤー設定",      p_set_layer,      3, 3},
    {"位置固定",          p_freeze_pos,     4, 4},
    {"回転固定",          p_freeze_rot,     4, 4},
    {"スリープ設定",      p_set_sleep,      2, 2},
    {"スリープ中",        p_is_sleeping,    1, 1},
    {"タグ設定",          p_set_tag,        2, 2},
    {"タグ取得",          p_get_tag,        1, 1},
    /* コライダー */
    {"球コライダー",      p_col_sphere,     2, 2},
    {"箱コライダー",      p_col_box,        4, 4},
    {"カプセルコライダー", p_col_capsule,    3, 3},
    {"コライダーオフセット",p_col_offset,   4, 4},
    /* 力 */
    {"力追加",            p_add_force,      4, 4},
    {"衝撃追加",          p_add_impulse,    4, 4},
    {"トルク追加",        p_add_torque,     4, 4},
    {"点力追加",          p_add_force_at,   7, 7},
    {"力消去",            p_clear_forces,   1, 1},
    /* ジョイント */
    {"固定ジョイント",    p_joint_fixed,    2, 2},
    {"距離ジョイント",    p_joint_distance, 4, 4},
    {"バネジョイント",    p_joint_spring,   5, 5},
    {"ヒンジジョイント",  p_joint_hinge,    7, 7},
    {"ジョイント削除",    p_joint_destroy,  1, 1},
    /* レイキャスト / クエリ */
    {"レイキャスト",      p_raycast,        6, 8},
    {"球キャスト",        p_sphere_cast,    7, 9},
    {"球重複検出",        p_overlap_sphere, 4, 5},
    {"箱重複検出",        p_overlap_box,    6, 7},
    {"衝突確認",          p_check_col,      2, 2},
    {"トリガー確認",      p_check_trig,     2, 2},
    /* ユーティリティ */
    {"前方向取得",        p_get_forward,    1, 1},
    {"上方向取得",        p_get_up,         1, 1},
    {"右方向取得",        p_get_right,      1, 1},
    {"速さ取得",          p_get_speed,      1, 1},
    {"接触数取得",        p_contact_count,  0, 0},
    {"接触剛体取得",      p_contact_bodies, 1, 1},
    {"接触情報取得",      p_contact_info,   1, 1},
};

HAJIMU_PLUGIN_EXPORT HajimuPluginInfo* hajimu_plugin_init(void){
    static HajimuPluginInfo info = {
        .name         = "jp_engine_physics",
        .version      = "1.0.0",
        .author       = "jp-physics contributors",
        .description  = "Unity レベル物理エンジン - 剛体/衝突/ジョイント/レイキャスト",
        .functions    = functions,
        .function_count = sizeof(functions)/sizeof(functions[0]),
    };
    return &info;
}
