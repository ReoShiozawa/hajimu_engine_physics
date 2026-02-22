/*  phys.h  ─  jp-physics  v1.0.0
 *  Unity レベル物理エンジン (はじむ言語プラグイン用)
 *  剛体・衝突・ジョイント・レイキャスト
 */
#pragma once
#include <stdint.h>

/* ─── 定数 ──────────────────────────────────────────── */
#define PHYS_MAX_BODIES    512
#define PHYS_MAX_CONTACTS  2048
#define PHYS_MAX_JOINTS    256
#define PHYS_MAX_LAYERS    32
#define PHYS_MAX_OVERLAP   64
#define PHYS_SLEEP_THRESH  0.02f    /* スリープ速度閾値 */
#define PHYS_SLEEP_FRAMES  60       /* スリープ判定フレーム数 */

/* ─── コライダー形状 ────────────────────────────────── */
typedef enum {
    PHYS_COL_NONE    = 0,
    PHYS_COL_SPHERE  = 1,
    PHYS_COL_BOX     = 2,
    PHYS_COL_CAPSULE = 3,
} PhysColType;

/* ─── 剛体タイプ ────────────────────────────────────── */
typedef enum {
    PHYS_BODY_DYNAMIC    = 0,
    PHYS_BODY_KINEMATIC  = 1,
    PHYS_BODY_STATIC     = 2,
} PhysBodyType;

/* ─── ジョイントタイプ ──────────────────────────────── */
typedef enum {
    PHYS_JOINT_FIXED    = 0,
    PHYS_JOINT_DISTANCE = 1,
    PHYS_JOINT_SPRING   = 2,
    PHYS_JOINT_HINGE    = 3,
} PhysJointType;

/* ─── ベクトル / Quaternion ─────────────────────────── */
typedef struct { float x, y, z;       } PhysVec3;
typedef struct { float x, y, z, w;    } PhysQuat;
typedef struct { float m[16];         } PhysMat4;

/* ─── 剛体 ──────────────────────────────────────────── */
typedef struct {
    /* 位置 / 姿勢 */
    PhysVec3    pos;
    PhysQuat    rot;
    /* 線形 / 角速度 */
    PhysVec3    vel;
    PhysVec3    omega;
    /* 累積力 / トルク */
    PhysVec3    force;
    PhysVec3    torque;
    /* 質量特性 */
    float       mass;
    float       inv_mass;
    PhysVec3    inv_inertia;   /* ローカル慣性テンソル(逆数) */
    /* 物理マテリアル */
    float       restitution;   /* 反発係数  0〜1 */
    float       friction;      /* 動摩擦係数 */
    float       static_friction;
    /* 減衰 */
    float       lin_drag;
    float       ang_drag;
    /* コライダー */
    PhysColType col_type;
    PhysVec3    col_size;      /* sphere:(r,0,0) box:(hx,hy,hz) capsule:(r,h,0) */
    PhysVec3    col_offset;    /* ローカルオフセット */
    /* タイプ / 状態 */
    PhysBodyType body_type;
    int         is_trigger;
    int         is_sleeping;
    int         sleep_counter;
    /* レイヤー */
    int         layer;         /* 0〜31 */
    uint32_t    layer_mask;    /* 衝突レイヤービットマスク */
    /* 軸固定 (Unity: Constraints) */
    int         freeze_pos;    /* bit0=X bit1=Y bit2=Z */
    int         freeze_rot;    /* bit0=X bit1=Y bit2=Z */
    /* ユーザーデータ */
    int         user_tag;
    /* 内部 */
    PhysVec3    aabb_min, aabb_max;  /* ワールドAABB */
    int         active;
} PhysBody;

/* ─── 衝突情報 ──────────────────────────────────────── */
typedef struct {
    int       a, b;           /* ボディインデックス */
    PhysVec3  normal;         /* a→b 法線 */
    float     depth;          /* 貫通深度 */
    PhysVec3  contact[2];     /* 接触点 */
    int       n_contacts;
    int       is_trigger;
} PhysContact;

/* ─── ジョイント ────────────────────────────────────── */
typedef struct {
    PhysJointType type;
    int       a, b;
    /* アンカー (ローカル空間) */
    PhysVec3  anchor_a;
    PhysVec3  anchor_b;
    /* 距離ジョイント */
    float     min_dist;
    float     max_dist;
    /* バネジョイント */
    float     rest_length;
    float     spring_k;
    float     spring_damper;
    /* ヒンジジョイント */
    PhysVec3  axis;           /* ローカル回転軸(a空間) */
    float     min_angle;      /* 度 */
    float     max_angle;      /* 度 */
    int       active;
} PhysJoint;

/* ─── レイヒット ────────────────────────────────────── */
typedef struct {
    int       hit;
    int       body_id;
    float     t;
    PhysVec3  point;
    PhysVec3  normal;
} PhysRayHit;

/* ─── ワールド ──────────────────────────────────────── */
typedef struct {
    PhysBody      bodies[PHYS_MAX_BODIES];
    PhysJoint     joints[PHYS_MAX_JOINTS];
    PhysContact   contacts[PHYS_MAX_CONTACTS];
    int           n_contacts;
    PhysVec3      gravity;
    float         fixed_dt;          /* サブステップ時間 */
    int           substeps;
    float         sleep_thresh;
    /* レイヤー衝突マトリクス (32×32 bit) */
    uint32_t      layer_matrix[PHYS_MAX_LAYERS];
    /* 衝突コールバックバッファ */
    int           trigger_pairs[PHYS_MAX_CONTACTS][2];
    int           n_triggers;
} PhysWorld;

/* ─── API ───────────────────────────────────────────── */
#ifdef __cplusplus
extern "C" {
#endif

PhysWorld* phys_world_create(void);
void       phys_world_destroy(PhysWorld* w);
void       phys_world_step(PhysWorld* w, float dt);
void       phys_set_gravity(PhysWorld* w, float x, float y, float z);
void       phys_set_substeps(PhysWorld* w, int n);
void       phys_set_layer_collision(PhysWorld* w, int la, int lb, int enabled);

/* 剛体 */
int        phys_body_create(PhysWorld* w);
void       phys_body_destroy(PhysWorld* w, int id);
void       phys_body_set_pos(PhysWorld* w, int id, float x, float y, float z);
void       phys_body_get_pos(PhysWorld* w, int id, float* x, float* y, float* z);
void       phys_body_set_rot_euler(PhysWorld* w, int id, float rx, float ry, float rz);
void       phys_body_get_rot_euler(PhysWorld* w, int id, float* rx, float* ry, float* rz);
void       phys_body_set_vel(PhysWorld* w, int id, float x, float y, float z);
void       phys_body_get_vel(PhysWorld* w, int id, float* x, float* y, float* z);
void       phys_body_set_omega(PhysWorld* w, int id, float x, float y, float z);
void       phys_body_get_omega(PhysWorld* w, int id, float* x, float* y, float* z);
void       phys_body_set_mass(PhysWorld* w, int id, float mass);
void       phys_body_set_drag(PhysWorld* w, int id, float lin, float ang);
void       phys_body_set_material(PhysWorld* w, int id, float restitution, float friction, float static_fric);
void       phys_body_set_type(PhysWorld* w, int id, PhysBodyType t);
void       phys_body_set_trigger(PhysWorld* w, int id, int enabled);
void       phys_body_set_layer(PhysWorld* w, int id, int layer, uint32_t mask);
void       phys_body_freeze_pos(PhysWorld* w, int id, int x, int y, int z);
void       phys_body_freeze_rot(PhysWorld* w, int id, int x, int y, int z);
void       phys_body_set_sleep(PhysWorld* w, int id, int sleeping);
int        phys_body_is_sleeping(PhysWorld* w, int id);
void       phys_body_set_tag(PhysWorld* w, int id, int tag);
int        phys_body_get_tag(PhysWorld* w, int id);

/* コライダー */
void       phys_col_sphere(PhysWorld* w, int id, float radius);
void       phys_col_box(PhysWorld* w, int id, float hx, float hy, float hz);
void       phys_col_capsule(PhysWorld* w, int id, float radius, float height);
void       phys_col_offset(PhysWorld* w, int id, float ox, float oy, float oz);

/* 力 */
void       phys_add_force(PhysWorld* w, int id, float fx, float fy, float fz);
void       phys_add_impulse(PhysWorld* w, int id, float ix, float iy, float iz);
void       phys_add_torque(PhysWorld* w, int id, float tx, float ty, float tz);
void       phys_add_force_at(PhysWorld* w, int id, float fx, float fy, float fz, float px, float py, float pz);
void       phys_clear_forces(PhysWorld* w, int id);

/* ジョイント */
int        phys_joint_fixed(PhysWorld* w, int a, int b);
int        phys_joint_distance(PhysWorld* w, int a, int b, float min_d, float max_d);
int        phys_joint_spring(PhysWorld* w, int a, int b, float rest, float k, float damper);
int        phys_joint_hinge(PhysWorld* w, int a, int b, float ax, float ay, float az, float min_deg, float max_deg);
void       phys_joint_destroy(PhysWorld* w, int jid);

/* レイキャスト / クエリ */
PhysRayHit phys_raycast(PhysWorld* w, float ox, float oy, float oz, float dx, float dy, float dz, float max_dist, uint32_t layer_mask);
PhysRayHit phys_sphere_cast(PhysWorld* w, float ox, float oy, float oz, float dx, float dy, float dz, float radius, float max_dist, uint32_t layer_mask);
int        phys_overlap_sphere(PhysWorld* w, float cx, float cy, float cz, float radius, uint32_t layer_mask, int* out, int max_out);
int        phys_overlap_box(PhysWorld* w, float cx, float cy, float cz, float hx, float hy, float hz, uint32_t layer_mask, int* out, int max_out);
int        phys_check_collision(PhysWorld* w, int a, int b);
int        phys_check_trigger(PhysWorld* w, int a, int b);

/* ユーティリティ */
void       phys_body_get_forward(PhysWorld* w, int id, float* x, float* y, float* z);
void       phys_body_get_up(PhysWorld* w, int id, float* x, float* y, float* z);
void       phys_body_get_right(PhysWorld* w, int id, float* x, float* y, float* z);
float      phys_body_get_speed(PhysWorld* w, int id);
int        phys_get_contact_count(PhysWorld* w);
int        phys_get_contact_bodies(PhysWorld* w, int ci, int* a, int* b);
void       phys_get_contact_info(PhysWorld* w, int ci, float* nx, float* ny, float* nz, float* depth);

#ifdef __cplusplus
}
#endif
