/*  phys.c  Part 1 — jp-physics v1.0.0
 *  数学ユーティリティ / ワールド / 剛体 / コライダー管理
 */
#include "phys.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdio.h>

/* ═══════════════════════════════════════════════════════
 *  数学ユーティリティ
 * ═══════════════════════════════════════════════════════ */
#define P_PI  3.14159265358979323846f
#define P_DEG2RAD (P_PI / 180.0f)
#define P_RAD2DEG (180.0f / P_PI)
#define P_EPS  1e-6f

static inline float pv3_dot(PhysVec3 a, PhysVec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static inline PhysVec3 pv3_cross(PhysVec3 a, PhysVec3 b) {
    return (PhysVec3){ a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x };
}
static inline PhysVec3 pv3_add(PhysVec3 a, PhysVec3 b) {
    return (PhysVec3){ a.x+b.x, a.y+b.y, a.z+b.z };
}
static inline PhysVec3 pv3_sub(PhysVec3 a, PhysVec3 b) {
    return (PhysVec3){ a.x-b.x, a.y-b.y, a.z-b.z };
}
static inline PhysVec3 pv3_scale(PhysVec3 a, float s) {
    return (PhysVec3){ a.x*s, a.y*s, a.z*s };
}
static inline PhysVec3 pv3_neg(PhysVec3 a) {
    return (PhysVec3){ -a.x, -a.y, -a.z };
}
static inline float pv3_len(PhysVec3 a) {
    return sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);
}
static inline float pv3_len2(PhysVec3 a) {
    return a.x*a.x + a.y*a.y + a.z*a.z;
}
static inline PhysVec3 pv3_norm(PhysVec3 a) {
    float l = pv3_len(a);
    if (l < P_EPS) return (PhysVec3){0,0,0};
    return pv3_scale(a, 1.0f/l);
}
static inline PhysVec3 pv3_lerp(PhysVec3 a, PhysVec3 b, float t) {
    return (PhysVec3){ a.x+(b.x-a.x)*t, a.y+(b.y-a.y)*t, a.z+(b.z-a.z)*t };
}
static inline PhysVec3 pv3_min3(PhysVec3 a, PhysVec3 b) {
    return (PhysVec3){ a.x<b.x?a.x:b.x, a.y<b.y?a.y:b.y, a.z<b.z?a.z:b.z };
}
static inline PhysVec3 pv3_max3(PhysVec3 a, PhysVec3 b) {
    return (PhysVec3){ a.x>b.x?a.x:b.x, a.y>b.y?a.y:b.y, a.z>b.z?a.z:b.z };
}
static inline PhysVec3 pv3_clamp(PhysVec3 v, PhysVec3 lo, PhysVec3 hi) {
    return pv3_max3(lo, pv3_min3(v, hi));
}
static inline float fclampf(float v, float lo, float hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

/* ─── Quaternion ─────────────────────────────────────── */
static inline PhysQuat pq_identity(void) {
    return (PhysQuat){0,0,0,1};
}
static inline PhysQuat pq_mul(PhysQuat a, PhysQuat b) {
    return (PhysQuat){
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
    };
}
static inline PhysQuat pq_conj(PhysQuat q) {
    return (PhysQuat){-q.x,-q.y,-q.z,q.w};
}
static inline PhysQuat pq_norm(PhysQuat q) {
    float l = sqrtf(q.x*q.x+q.y*q.y+q.z*q.z+q.w*q.w);
    if (l < P_EPS) return pq_identity();
    float inv = 1.0f/l;
    return (PhysQuat){q.x*inv,q.y*inv,q.z*inv,q.w*inv};
}
static inline PhysVec3 pq_rot(PhysQuat q, PhysVec3 v) {
    /* p' = q * p * q^{-1} */
    PhysQuat p = {v.x, v.y, v.z, 0};
    PhysQuat r = pq_mul(pq_mul(q, p), pq_conj(q));
    return (PhysVec3){r.x, r.y, r.z};
}
/* オイラー角 (度) → クォータニオン (YXZ順) */
static PhysQuat pq_from_euler(float rx, float ry, float rz) {
    float hr = rx * P_DEG2RAD * 0.5f;
    float hp = ry * P_DEG2RAD * 0.5f;
    float hyw = rz * P_DEG2RAD * 0.5f;
    float cr = cosf(hr), sr = sinf(hr);
    float cp = cosf(hp), sp = sinf(hp);
    float cy = cosf(hyw), sy = sinf(hyw);
    return (PhysQuat){
        sr*cp*cy - cr*sp*sy,
        cr*sp*cy + sr*cp*sy,
        cr*cp*sy - sr*sp*cy,
        cr*cp*cy + sr*sp*sy
    };
}
/* クォータニオン → オイラー角 (度,YXZ) */
static void pq_to_euler(PhysQuat q, float* rx, float* ry, float* rz) {
    float sinr_cosp = 2.0f*(q.w*q.x + q.y*q.z);
    float cosr_cosp = 1.0f - 2.0f*(q.x*q.x + q.y*q.y);
    *rx = atan2f(sinr_cosp, cosr_cosp) * P_RAD2DEG;
    float sinp = 2.0f*(q.w*q.y - q.z*q.x);
    if (fabsf(sinp) >= 1.0f) *ry = copysignf(90.0f, sinp);
    else *ry = asinf(sinp) * P_RAD2DEG;
    float siny_cosp = 2.0f*(q.w*q.z + q.x*q.y);
    float cosy_cosp = 1.0f - 2.0f*(q.y*q.y + q.z*q.z);
    *rz = atan2f(siny_cosp, cosy_cosp) * P_RAD2DEG;
}
/* 角速度を使ってクォータニオンを積分 */
static PhysQuat pq_integrate(PhysQuat q, PhysVec3 omega, float dt) {
    float angle = pv3_len(omega) * dt;
    if (angle < P_EPS) return q;
    PhysVec3 axis = pv3_norm(omega);
    float ha = angle * 0.5f;
    PhysQuat dq = {
        axis.x * sinf(ha),
        axis.y * sinf(ha),
        axis.z * sinf(ha),
        cosf(ha)
    };
    return pq_norm(pq_mul(dq, q));
}

/* ─── 慣性テンソル (球/箱/カプセル ローカル逆数) ───────── */
static PhysVec3 compute_inv_inertia(PhysColType t, PhysVec3 sz, float mass) {
    if (mass <= 0.0f) return (PhysVec3){0,0,0};
    float inv_m = 1.0f / mass;
    switch (t) {
    case PHYS_COL_SPHERE: {
        float r = sz.x;
        float I = (2.0f/5.0f) * mass * r * r;
        float inv_I = 1.0f / I;
        return (PhysVec3){inv_I, inv_I, inv_I};
    }
    case PHYS_COL_BOX: {
        float hx=sz.x, hy=sz.y, hz=sz.z;
        float Ix = (1.0f/12.0f)*mass*(4.0f*hy*hy + 4.0f*hz*hz);
        float Iy = (1.0f/12.0f)*mass*(4.0f*hx*hx + 4.0f*hz*hz);
        float Iz = (1.0f/12.0f)*mass*(4.0f*hx*hx + 4.0f*hy*hy);
        return (PhysVec3){1.0f/Ix, 1.0f/Iy, 1.0f/Iz};
    }
    case PHYS_COL_CAPSULE: {
        float r = sz.x, h = sz.y;
        float ms = mass * (r*r*4.0f/3.0f*P_PI) / (r*r*4.0f/3.0f*P_PI + r*r*P_PI*h);
        float mc = mass - ms;
        float Ixx = mc*(h*h/12.0f + r*r/4.0f) + ms*(2.0f*r*r/5.0f + h*h/2.0f + 3.0f*h*r/8.0f);
        float Iyy = (mc * r*r / 2.0f) + ms * (2.0f*r*r/5.0f);
        (void)inv_m;
        return (PhysVec3){1.0f/Ixx, 1.0f/Iyy, 1.0f/Ixx};
    }
    default:
        return (PhysVec3){inv_m, inv_m, inv_m};
    }
}

/* ─── AABB 計算 ──────────────────────────────────────── */
static void compute_aabb(PhysBody* b) {
    PhysVec3 center = pv3_add(b->pos, pq_rot(b->rot, b->col_offset));
    float pad = 0.0f;
    switch (b->col_type) {
    case PHYS_COL_SPHERE: {
        float r = b->col_size.x + pad;
        b->aabb_min = (PhysVec3){center.x-r, center.y-r, center.z-r};
        b->aabb_max = (PhysVec3){center.x+r, center.y+r, center.z+r};
        break;
    }
    case PHYS_COL_BOX: {
        /* 回転したAABBの外接AABB */
        float hx=b->col_size.x, hy=b->col_size.y, hz=b->col_size.z;
        PhysVec3 ax = pq_rot(b->rot, (PhysVec3){hx,0,0});
        PhysVec3 ay = pq_rot(b->rot, (PhysVec3){0,hy,0});
        PhysVec3 az = pq_rot(b->rot, (PhysVec3){0,0,hz});
        float ex = fabsf(ax.x)+fabsf(ay.x)+fabsf(az.x);
        float ey = fabsf(ax.y)+fabsf(ay.y)+fabsf(az.y);
        float ez = fabsf(ax.z)+fabsf(ay.z)+fabsf(az.z);
        b->aabb_min = (PhysVec3){center.x-ex-pad, center.y-ey-pad, center.z-ez-pad};
        b->aabb_max = (PhysVec3){center.x+ex+pad, center.y+ey+pad, center.z+ez+pad};
        break;
    }
    case PHYS_COL_CAPSULE: {
        float r = b->col_size.x, h = b->col_size.y;
        /* カプセルはY軸方向 */
        PhysVec3 up = pq_rot(b->rot, (PhysVec3){0,1,0});
        PhysVec3 p0 = pv3_add(center, pv3_scale(up,  h*0.5f));
        PhysVec3 p1 = pv3_add(center, pv3_scale(up, -h*0.5f));
        b->aabb_min = pv3_sub(pv3_min3(p0,p1), (PhysVec3){r,r,r});
        b->aabb_max = pv3_add(pv3_max3(p0,p1), (PhysVec3){r,r,r});
        break;
    }
    default:
        b->aabb_min = pv3_sub(center, (PhysVec3){0.5f,0.5f,0.5f});
        b->aabb_max = pv3_add(center, (PhysVec3){0.5f,0.5f,0.5f});
        break;
    }
}

static int aabb_overlap(PhysBody* a, PhysBody* b) {
    return (a->aabb_min.x <= b->aabb_max.x && a->aabb_max.x >= b->aabb_min.x &&
            a->aabb_min.y <= b->aabb_max.y && a->aabb_max.y >= b->aabb_min.y &&
            a->aabb_min.z <= b->aabb_max.z && a->aabb_max.z >= b->aabb_min.z);
}

/* ═══════════════════════════════════════════════════════
 *  ワールド管理
 * ═══════════════════════════════════════════════════════ */
PhysWorld* phys_world_create(void) {
    PhysWorld* w = (PhysWorld*)calloc(1, sizeof(PhysWorld));
    if (!w) return NULL;
    w->gravity    = (PhysVec3){0, -9.81f, 0};
    w->fixed_dt   = 1.0f / 60.0f;
    w->substeps   = 4;
    w->sleep_thresh = PHYS_SLEEP_THRESH;
    /* デフォルト: 全レイヤー同士は衝突 */
    for (int i = 0; i < PHYS_MAX_LAYERS; i++)
        w->layer_matrix[i] = 0xFFFFFFFF;
    return w;
}
void phys_world_destroy(PhysWorld* w) {
    if (w) free(w);
}
void phys_set_gravity(PhysWorld* w, float x, float y, float z) {
    if (!w) return;
    w->gravity = (PhysVec3){x,y,z};
}
void phys_set_substeps(PhysWorld* w, int n) {
    if (!w || n < 1) return;
    w->substeps = n > 16 ? 16 : n;
}
void phys_set_layer_collision(PhysWorld* w, int la, int lb, int enabled) {
    if (!w || la<0||la>=PHYS_MAX_LAYERS||lb<0||lb>=PHYS_MAX_LAYERS) return;
    if (enabled) {
        w->layer_matrix[la] |=  (1u << lb);
        w->layer_matrix[lb] |=  (1u << la);
    } else {
        w->layer_matrix[la] &= ~(1u << lb);
        w->layer_matrix[lb] &= ~(1u << la);
    }
}

/* ═══════════════════════════════════════════════════════
 *  剛体管理
 * ═══════════════════════════════════════════════════════ */
int phys_body_create(PhysWorld* w) {
    if (!w) return -1;
    for (int i = 0; i < PHYS_MAX_BODIES; i++) {
        if (!w->bodies[i].active) {
            PhysBody* b = &w->bodies[i];
            memset(b, 0, sizeof(PhysBody));
            b->active      = 1;
            b->rot         = pq_identity();
            b->mass        = 1.0f;
            b->inv_mass    = 1.0f;
            b->restitution = 0.3f;
            b->friction    = 0.5f;
            b->static_friction = 0.6f;
            b->lin_drag    = 0.05f;
            b->ang_drag    = 0.08f;
            b->body_type   = PHYS_BODY_DYNAMIC;
            b->layer       = 0;
            b->layer_mask  = 0xFFFFFFFF;
            b->col_type    = PHYS_COL_SPHERE;
            b->col_size    = (PhysVec3){0.5f, 0.5f, 0.5f};
            b->inv_inertia = compute_inv_inertia(PHYS_COL_SPHERE, b->col_size, 1.0f);
            return i;
        }
    }
    return -1;
}
void phys_body_destroy(PhysWorld* w, int id) {
    if (!w||id<0||id>=PHYS_MAX_BODIES) return;
    w->bodies[id].active = 0;
    /* 関係ジョイントも無効化 */
    for (int j = 0; j < PHYS_MAX_JOINTS; j++) {
        if (w->joints[j].active && (w->joints[j].a==id||w->joints[j].b==id))
            w->joints[j].active = 0;
    }
}

/* 位置 / 回転 */
void phys_body_set_pos(PhysWorld* w, int id, float x, float y, float z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].pos = (PhysVec3){x,y,z};
}
void phys_body_get_pos(PhysWorld* w, int id, float* x, float* y, float* z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) { *x=*y=*z=0; return; }
    *x = w->bodies[id].pos.x;
    *y = w->bodies[id].pos.y;
    *z = w->bodies[id].pos.z;
}
void phys_body_set_rot_euler(PhysWorld* w, int id, float rx, float ry, float rz) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].rot = pq_from_euler(rx, ry, rz);
}
void phys_body_get_rot_euler(PhysWorld* w, int id, float* rx, float* ry, float* rz) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) { *rx=*ry=*rz=0; return; }
    pq_to_euler(w->bodies[id].rot, rx, ry, rz);
}

/* 速度 */
void phys_body_set_vel(PhysWorld* w, int id, float x, float y, float z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].vel = (PhysVec3){x,y,z};
}
void phys_body_get_vel(PhysWorld* w, int id, float* x, float* y, float* z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) { *x=*y=*z=0; return; }
    *x=w->bodies[id].vel.x; *y=w->bodies[id].vel.y; *z=w->bodies[id].vel.z;
}
void phys_body_set_omega(PhysWorld* w, int id, float x, float y, float z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].omega = (PhysVec3){x,y,z};
}
void phys_body_get_omega(PhysWorld* w, int id, float* x, float* y, float* z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) { *x=*y=*z=0; return; }
    *x=w->bodies[id].omega.x; *y=w->bodies[id].omega.y; *z=w->bodies[id].omega.z;
}

/* 質量 */
void phys_body_set_mass(PhysWorld* w, int id, float mass) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    PhysBody* b = &w->bodies[id];
    b->mass = mass;
    b->inv_mass = (mass > 0.0f) ? 1.0f/mass : 0.0f;
    b->inv_inertia = compute_inv_inertia(b->col_type, b->col_size, mass);
}
void phys_body_set_drag(PhysWorld* w, int id, float lin, float ang) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].lin_drag = lin;
    w->bodies[id].ang_drag = ang;
}
void phys_body_set_material(PhysWorld* w, int id, float res, float fric, float sfric) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].restitution = res;
    w->bodies[id].friction    = fric;
    w->bodies[id].static_friction = sfric;
}
void phys_body_set_type(PhysWorld* w, int id, PhysBodyType t) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    PhysBody* b = &w->bodies[id];
    b->body_type = t;
    if (t == PHYS_BODY_STATIC) { b->inv_mass = 0.0f; b->inv_inertia=(PhysVec3){0,0,0}; }
    else if (t == PHYS_BODY_KINEMATIC) { b->inv_mass = 0.0f; b->inv_inertia=(PhysVec3){0,0,0}; }
    else { phys_body_set_mass(w, id, b->mass); }
}
void phys_body_set_trigger(PhysWorld* w, int id, int en) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].is_trigger = en;
}
void phys_body_set_layer(PhysWorld* w, int id, int layer, uint32_t mask) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].layer = (layer >= 0 && layer < PHYS_MAX_LAYERS) ? layer : 0;
    w->bodies[id].layer_mask = mask;
}
void phys_body_freeze_pos(PhysWorld* w, int id, int x, int y, int z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].freeze_pos = (x?1:0)|(y?2:0)|(z?4:0);
}
void phys_body_freeze_rot(PhysWorld* w, int id, int x, int y, int z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].freeze_rot = (x?1:0)|(y?2:0)|(z?4:0);
}
void phys_body_set_sleep(PhysWorld* w, int id, int sleeping) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].is_sleeping  = sleeping;
    w->bodies[id].sleep_counter= sleeping ? PHYS_SLEEP_FRAMES : 0;
}
int phys_body_is_sleeping(PhysWorld* w, int id) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return 0;
    return w->bodies[id].is_sleeping;
}
void phys_body_set_tag(PhysWorld* w, int id, int tag) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].user_tag = tag;
}
int phys_body_get_tag(PhysWorld* w, int id) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return -1;
    return w->bodies[id].user_tag;
}

/* コライダー */
void phys_col_sphere(PhysWorld* w, int id, float r) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    PhysBody* b = &w->bodies[id];
    b->col_type = PHYS_COL_SPHERE;
    b->col_size = (PhysVec3){r, r, r};
    b->inv_inertia = compute_inv_inertia(PHYS_COL_SPHERE, b->col_size, b->mass);
}
void phys_col_box(PhysWorld* w, int id, float hx, float hy, float hz) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    PhysBody* b = &w->bodies[id];
    b->col_type = PHYS_COL_BOX;
    b->col_size = (PhysVec3){hx, hy, hz};
    b->inv_inertia = compute_inv_inertia(PHYS_COL_BOX, b->col_size, b->mass);
}
void phys_col_capsule(PhysWorld* w, int id, float r, float h) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    PhysBody* b = &w->bodies[id];
    b->col_type = PHYS_COL_CAPSULE;
    b->col_size = (PhysVec3){r, h, 0};
    b->inv_inertia = compute_inv_inertia(PHYS_COL_CAPSULE, b->col_size, b->mass);
}
void phys_col_offset(PhysWorld* w, int id, float ox, float oy, float oz) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].col_offset = (PhysVec3){ox, oy, oz};
}

/* 力 */
void phys_add_force(PhysWorld* w, int id, float fx, float fy, float fz) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    if (w->bodies[id].body_type != PHYS_BODY_DYNAMIC) return;
    w->bodies[id].force.x += fx;
    w->bodies[id].force.y += fy;
    w->bodies[id].force.z += fz;
}
void phys_add_impulse(PhysWorld* w, int id, float ix, float iy, float iz) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    PhysBody* b = &w->bodies[id];
    if (b->body_type != PHYS_BODY_DYNAMIC) return;
    b->vel.x += ix * b->inv_mass;
    b->vel.y += iy * b->inv_mass;
    b->vel.z += iz * b->inv_mass;
    b->is_sleeping = 0;
    b->sleep_counter = 0;
}
void phys_add_torque(PhysWorld* w, int id, float tx, float ty, float tz) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    if (w->bodies[id].body_type != PHYS_BODY_DYNAMIC) return;
    w->bodies[id].torque.x += tx;
    w->bodies[id].torque.y += ty;
    w->bodies[id].torque.z += tz;
}
void phys_add_force_at(PhysWorld* w, int id, float fx, float fy, float fz, float px, float py, float pz) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    PhysBody* b = &w->bodies[id];
    if (b->body_type != PHYS_BODY_DYNAMIC) return;
    b->force.x += fx; b->force.y += fy; b->force.z += fz;
    PhysVec3 r = pv3_sub((PhysVec3){px,py,pz}, b->pos);
    PhysVec3 f = (PhysVec3){fx,fy,fz};
    PhysVec3 tau = pv3_cross(r, f);
    b->torque.x += tau.x; b->torque.y += tau.y; b->torque.z += tau.z;
}
void phys_clear_forces(PhysWorld* w, int id) {
    if (!w||id<0||id>=PHYS_MAX_BODIES||!w->bodies[id].active) return;
    w->bodies[id].force  = (PhysVec3){0,0,0};
    w->bodies[id].torque = (PhysVec3){0,0,0};
}

/* ユーティリティ */
void phys_body_get_forward(PhysWorld* w, int id, float* x, float* y, float* z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES) { *x=0;*y=0;*z=-1; return; }
    PhysVec3 f = pq_rot(w->bodies[id].rot, (PhysVec3){0,0,-1});
    *x=f.x; *y=f.y; *z=f.z;
}
void phys_body_get_up(PhysWorld* w, int id, float* x, float* y, float* z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES) { *x=0;*y=1;*z=0; return; }
    PhysVec3 u = pq_rot(w->bodies[id].rot, (PhysVec3){0,1,0});
    *x=u.x; *y=u.y; *z=u.z;
}
void phys_body_get_right(PhysWorld* w, int id, float* x, float* y, float* z) {
    if (!w||id<0||id>=PHYS_MAX_BODIES) { *x=1;*y=0;*z=0; return; }
    PhysVec3 r = pq_rot(w->bodies[id].rot, (PhysVec3){1,0,0});
    *x=r.x; *y=r.y; *z=r.z;
}
float phys_body_get_speed(PhysWorld* w, int id) {
    if (!w||id<0||id>=PHYS_MAX_BODIES) return 0.0f;
    return pv3_len(w->bodies[id].vel);
}
int phys_get_contact_count(PhysWorld* w) {
    if (!w) return 0;
    return w->n_contacts;
}
int phys_get_contact_bodies(PhysWorld* w, int ci, int* a, int* b) {
    if (!w||ci<0||ci>=w->n_contacts) return 0;
    *a = w->contacts[ci].a;
    *b = w->contacts[ci].b;
    return 1;
}
void phys_get_contact_info(PhysWorld* w, int ci, float* nx, float* ny, float* nz, float* depth) {
    if (!w||ci<0||ci>=w->n_contacts) { *nx=0;*ny=1;*nz=0;*depth=0; return; }
    *nx=w->contacts[ci].normal.x;
    *ny=w->contacts[ci].normal.y;
    *nz=w->contacts[ci].normal.z;
    *depth=w->contacts[ci].depth;
}

/* ─── ジョイント作成 ─────────────────────────────────── */
static int alloc_joint(PhysWorld* w) {
    for (int i = 0; i < PHYS_MAX_JOINTS; i++)
        if (!w->joints[i].active) return i;
    return -1;
}
int phys_joint_fixed(PhysWorld* w, int a, int b) {
    int ji = alloc_joint(w); if (ji<0) return -1;
    PhysJoint* j = &w->joints[ji];
    memset(j, 0, sizeof(PhysJoint));
    j->type=PHYS_JOINT_FIXED; j->a=a; j->b=b; j->active=1;
    return ji;
}
int phys_joint_distance(PhysWorld* w, int a, int b, float mn, float mx) {
    int ji = alloc_joint(w); if (ji<0) return -1;
    PhysJoint* j = &w->joints[ji];
    memset(j, 0, sizeof(PhysJoint));
    j->type=PHYS_JOINT_DISTANCE; j->a=a; j->b=b;
    j->min_dist=mn; j->max_dist=mx; j->active=1;
    return ji;
}
int phys_joint_spring(PhysWorld* w, int a, int b, float rest, float k, float damper) {
    int ji = alloc_joint(w); if (ji<0) return -1;
    PhysJoint* j = &w->joints[ji];
    memset(j, 0, sizeof(PhysJoint));
    j->type=PHYS_JOINT_SPRING; j->a=a; j->b=b;
    j->rest_length=rest; j->spring_k=k; j->spring_damper=damper; j->active=1;
    return ji;
}
int phys_joint_hinge(PhysWorld* w, int a, int b, float ax, float ay, float az, float mn, float mx) {
    int ji = alloc_joint(w); if (ji<0) return -1;
    PhysJoint* j = &w->joints[ji];
    memset(j, 0, sizeof(PhysJoint));
    j->type=PHYS_JOINT_HINGE; j->a=a; j->b=b;
    j->axis=pv3_norm((PhysVec3){ax,ay,az});
    j->min_angle=mn; j->max_angle=mx; j->active=1;
    return ji;
}
void phys_joint_destroy(PhysWorld* w, int ji) {
    if (!w||ji<0||ji>=PHYS_MAX_JOINTS) return;
    w->joints[ji].active = 0;
}
/* ═══════════════════════════════════════════════════════
 *  phys.c  Part 2 — ナロウフェーズ衝突検出
 *  Sphere-Sphere / Sphere-Box / Box-Box /
 *  Sphere-Capsule / Capsule-Capsule / Capsule-Box
 * ═══════════════════════════════════════════════════════ */

/* ─── 最近傍点ユーティリティ ─────────────────────────── */

/* 点→線分の最近傍パラメータ t (0〜1) */
static float seg_closest_t(PhysVec3 p, PhysVec3 a, PhysVec3 b) {
    PhysVec3 ab = pv3_sub(b,a);
    float len2 = pv3_len2(ab);
    if (len2 < P_EPS) return 0.0f;
    return fclampf(pv3_dot(pv3_sub(p,a), ab)/len2, 0.0f, 1.0f);
}
/* 2線分間の最近傍点ペア */
static void seg_seg_closest(PhysVec3 a0, PhysVec3 a1,
                             PhysVec3 b0, PhysVec3 b1,
                             PhysVec3* ca, PhysVec3* cb) {
    PhysVec3 da = pv3_sub(a1,a0), db = pv3_sub(b1,b0);
    PhysVec3 r  = pv3_sub(a0,b0);
    float a = pv3_dot(da,da), e = pv3_dot(db,db), f = pv3_dot(db,r);
    float s, t;
    if (a < P_EPS && e < P_EPS) { *ca=a0; *cb=b0; return; }
    if (a < P_EPS) { t=fclampf(f/e,0,1); s=0; }
    else {
        float c = pv3_dot(da,r);
        if (e < P_EPS) { s=fclampf(-c/a,0,1); t=0; }
        else {
            float b2 = pv3_dot(da,db), denom = a*e - b2*b2;
            if (fabsf(denom) > P_EPS) s=fclampf((b2*f-c*e)/denom,0,1);
            else s=0;
            t=fclampf((b2*s+f)/e,0,1);
            s=fclampf((b2*t-c)/a,0,1);
            t=fclampf((b2*s+f)/e,0,1);
        }
    }
    *ca = pv3_add(a0, pv3_scale(da,s));
    *cb = pv3_add(b0, pv3_scale(db,t));
}
/* 点→OBBの最近傍点 */
static PhysVec3 point_closest_on_obb(PhysVec3 p, PhysVec3 center,
                                      PhysVec3 ax, PhysVec3 ay, PhysVec3 az,
                                      float hx, float hy, float hz) {
    PhysVec3 d = pv3_sub(p, center);
    float tx = fclampf(pv3_dot(d,ax),-hx,hx);
    float ty = fclampf(pv3_dot(d,ay),-hy,hy);
    float tz = fclampf(pv3_dot(d,az),-hz,hz);
    return pv3_add(pv3_add(pv3_add(center, pv3_scale(ax,tx)),
                                            pv3_scale(ay,ty)),
                                            pv3_scale(az,tz));
}
/* カプセルの2端点を取得 */
static void capsule_endpoints(PhysBody* b, PhysVec3* p0, PhysVec3* p1) {
    PhysVec3 center = pv3_add(b->pos, pq_rot(b->rot, b->col_offset));
    PhysVec3 up = pq_rot(b->rot, (PhysVec3){0,1,0});
    float h2 = b->col_size.y * 0.5f;
    *p0 = pv3_add(center, pv3_scale(up,  h2));
    *p1 = pv3_add(center, pv3_scale(up, -h2));
}

/* ─── Sphere vs Sphere ───────────────────────────────── */
static int detect_sphere_sphere(PhysBody* a, PhysBody* b, PhysContact* c) {
    PhysVec3 ca = pv3_add(a->pos, pq_rot(a->rot, a->col_offset));
    PhysVec3 cb = pv3_add(b->pos, pq_rot(b->rot, b->col_offset));
    float ra = a->col_size.x, rb = b->col_size.x;
    PhysVec3 d = pv3_sub(cb, ca);
    float dist2 = pv3_len2(d);
    float r_sum = ra + rb;
    if (dist2 >= r_sum*r_sum) return 0;
    float dist = sqrtf(dist2);
    c->normal   = dist > P_EPS ? pv3_scale(d, 1.0f/dist) : (PhysVec3){0,1,0};
    c->depth    = r_sum - dist;
    c->contact[0] = pv3_add(ca, pv3_scale(c->normal, ra));
    c->n_contacts = 1;
    return 1;
}

/* ─── Sphere vs Box ──────────────────────────────────── */
static int detect_sphere_box(PhysBody* sa, PhysBody* ba_body, PhysContact* c) {
    PhysVec3 sc = pv3_add(sa->pos, pq_rot(sa->rot, sa->col_offset));
    float sr    = sa->col_size.x;
    PhysVec3 bc = pv3_add(ba_body->pos, pq_rot(ba_body->rot, ba_body->col_offset));
    float hx=ba_body->col_size.x, hy=ba_body->col_size.y, hz=ba_body->col_size.z;
    PhysVec3 bax = pq_rot(ba_body->rot,(PhysVec3){1,0,0});
    PhysVec3 bay = pq_rot(ba_body->rot,(PhysVec3){0,1,0});
    PhysVec3 baz = pq_rot(ba_body->rot,(PhysVec3){0,0,1});
    PhysVec3 closest = point_closest_on_obb(sc, bc, bax, bay, baz, hx, hy, hz);
    PhysVec3 d = pv3_sub(sc, closest);
    float dist2 = pv3_len2(d);
    if (dist2 >= sr*sr) return 0;
    float dist = sqrtf(dist2);
    c->normal   = dist > P_EPS ? pv3_scale(d, 1.0f/dist) : (PhysVec3){0,1,0};
    c->depth    = sr - dist;
    c->contact[0] = closest;
    c->n_contacts = 1;
    return 1;
}

/* ─── Box vs Box (SAT) ───────────────────────────────── */
static float box_overlap_on_axis(PhysBody* a, PhysBody* b, PhysVec3 axis) {
    /* 各OBBの半径(投影) */
    float ra =
        fabsf(pv3_dot(pq_rot(a->rot,(PhysVec3){1,0,0}), axis)) * a->col_size.x +
        fabsf(pv3_dot(pq_rot(a->rot,(PhysVec3){0,1,0}), axis)) * a->col_size.y +
        fabsf(pv3_dot(pq_rot(a->rot,(PhysVec3){0,0,1}), axis)) * a->col_size.z;
    float rb =
        fabsf(pv3_dot(pq_rot(b->rot,(PhysVec3){1,0,0}), axis)) * b->col_size.x +
        fabsf(pv3_dot(pq_rot(b->rot,(PhysVec3){0,1,0}), axis)) * b->col_size.y +
        fabsf(pv3_dot(pq_rot(b->rot,(PhysVec3){0,0,1}), axis)) * b->col_size.z;
    float ca_proj = pv3_dot(pv3_add(a->pos, pq_rot(a->rot, a->col_offset)), axis);
    float cb_proj = pv3_dot(pv3_add(b->pos, pq_rot(b->rot, b->col_offset)), axis);
    return (ra + rb) - fabsf(ca_proj - cb_proj);
}

static int detect_box_box(PhysBody* a, PhysBody* b, PhysContact* c) {
    /* 15軸 SAT (3+3+9) */
    PhysVec3 axes[15];
    int n = 0;
    PhysVec3 ax[3], bx[3];
    ax[0]=pq_rot(a->rot,(PhysVec3){1,0,0}); ax[1]=pq_rot(a->rot,(PhysVec3){0,1,0}); ax[2]=pq_rot(a->rot,(PhysVec3){0,0,1});
    bx[0]=pq_rot(b->rot,(PhysVec3){1,0,0}); bx[1]=pq_rot(b->rot,(PhysVec3){0,1,0}); bx[2]=pq_rot(b->rot,(PhysVec3){0,0,1});
    for (int i=0;i<3;i++) { axes[n++]=ax[i]; axes[n++]=bx[i]; }
    for (int i=0;i<3;i++)
        for (int j=0;j<3;j++) {
            PhysVec3 cr = pv3_cross(ax[i], bx[j]);
            if (pv3_len2(cr) > P_EPS) axes[n++] = pv3_norm(cr);
            else n++; /* padding — will be filtered by len check */
        }
    float min_overlap = FLT_MAX;
    PhysVec3 best_axis = {0,1,0};
    PhysVec3 d = pv3_sub(pv3_add(b->pos,pq_rot(b->rot,b->col_offset)),
                          pv3_add(a->pos,pq_rot(a->rot,a->col_offset)));
    for (int i = 0; i < n; i++) {
        if (pv3_len2(axes[i]) < P_EPS) continue;
        PhysVec3 ax_n = pv3_norm(axes[i]);
        float ov = box_overlap_on_axis(a, b, ax_n);
        if (ov <= 0.0f) return 0;  /* 分離軸発見 */
        if (ov < min_overlap) {
            min_overlap = ov;
            best_axis = pv3_dot(d, ax_n) < 0 ? pv3_neg(ax_n) : ax_n;
        }
    }
    c->normal   = best_axis;
    c->depth    = min_overlap;
    PhysVec3 ca = pv3_add(a->pos, pq_rot(a->rot, a->col_offset));
    c->contact[0] = pv3_add(ca, pv3_scale(best_axis, min_overlap*0.5f));
    c->n_contacts = 1;
    return 1;
}

/* ─── Sphere vs Capsule ──────────────────────────────── */
static int detect_sphere_capsule(PhysBody* sa, PhysBody* cap, PhysContact* c) {
    PhysVec3 sc = pv3_add(sa->pos, pq_rot(sa->rot, sa->col_offset));
    float sr    = sa->col_size.x;
    PhysVec3 cp0, cp1;
    capsule_endpoints(cap, &cp0, &cp1);
    float t = seg_closest_t(sc, cp0, cp1);
    PhysVec3 closest = pv3_lerp(cp0, cp1, t);
    PhysVec3 d = pv3_sub(sc, closest);
    float dist2 = pv3_len2(d);
    float r_sum = sr + cap->col_size.x;
    if (dist2 >= r_sum*r_sum) return 0;
    float dist = sqrtf(dist2);
    c->normal   = dist > P_EPS ? pv3_scale(d, 1.0f/dist) : (PhysVec3){0,1,0};
    c->depth    = r_sum - dist;
    c->contact[0] = pv3_add(closest, pv3_scale(c->normal, cap->col_size.x));
    c->n_contacts = 1;
    return 1;
}

/* ─── Capsule vs Capsule ─────────────────────────────── */
static int detect_capsule_capsule(PhysBody* a, PhysBody* b, PhysContact* c) {
    PhysVec3 a0,a1,b0,b1;
    capsule_endpoints(a, &a0, &a1);
    capsule_endpoints(b, &b0, &b1);
    PhysVec3 ca, cb;
    seg_seg_closest(a0,a1,b0,b1,&ca,&cb);
    PhysVec3 d = pv3_sub(ca, cb);
    float dist2 = pv3_len2(d);
    float r_sum = a->col_size.x + b->col_size.x;
    if (dist2 >= r_sum*r_sum) return 0;
    float dist = sqrtf(dist2);
    c->normal   = dist > P_EPS ? pv3_scale(d, 1.0f/dist) : (PhysVec3){0,1,0};
    c->depth    = r_sum - dist;
    c->contact[0] = pv3_add(cb, pv3_scale(c->normal, b->col_size.x));
    c->n_contacts = 1;
    return 1;
}

/* ─── Capsule vs Box ─────────────────────────────────── */
static int detect_capsule_box(PhysBody* cap, PhysBody* box, PhysContact* c) {
    PhysVec3 p0, p1;
    capsule_endpoints(cap, &p0, &p1);
    float r = cap->col_size.x;
    PhysVec3 bc = pv3_add(box->pos, pq_rot(box->rot, box->col_offset));
    float hx=box->col_size.x,hy=box->col_size.y,hz=box->col_size.z;
    PhysVec3 bax=pq_rot(box->rot,(PhysVec3){1,0,0});
    PhysVec3 bay=pq_rot(box->rot,(PhysVec3){0,1,0});
    PhysVec3 baz=pq_rot(box->rot,(PhysVec3){0,0,1});
    /* 線分をOBBローカル空間へ変換して最近傍点を探す */
    PhysVec3 best_cap={0}, best_box={0};
    float best_dist2 = FLT_MAX;
    int N = 8;
    for (int i = 0; i <= N; i++) {
        float t = (float)i / N;
        PhysVec3 pt = pv3_lerp(p0, p1, t);
        PhysVec3 cl = point_closest_on_obb(pt, bc, bax, bay, baz, hx, hy, hz);
        float d2 = pv3_len2(pv3_sub(pt, cl));
        if (d2 < best_dist2) {
            best_dist2 = d2; best_cap = pt; best_box = cl;
        }
    }
    if (best_dist2 >= r*r) return 0;
    float dist = sqrtf(best_dist2);
    PhysVec3 d = pv3_sub(best_cap, best_box);
    c->normal   = dist > P_EPS ? pv3_scale(d, 1.0f/dist) : (PhysVec3){0,1,0};
    c->depth    = r - dist;
    c->contact[0] = best_box;
    c->n_contacts = 1;
    return 1;
}

/* ─── ディスパッチ ────────────────────────────────────── */
static int detect_pair(PhysBody* a, PhysBody* b, PhysContact* c) {
    /* A が Sphere のケース */
    if (a->col_type==PHYS_COL_SPHERE && b->col_type==PHYS_COL_SPHERE)
        return detect_sphere_sphere(a, b, c);
    if (a->col_type==PHYS_COL_SPHERE && b->col_type==PHYS_COL_BOX)
        return detect_sphere_box(a, b, c);
    if (a->col_type==PHYS_COL_BOX    && b->col_type==PHYS_COL_SPHERE) {
        int r = detect_sphere_box(b, a, c);
        c->normal = pv3_neg(c->normal); return r;
    }
    if (a->col_type==PHYS_COL_BOX    && b->col_type==PHYS_COL_BOX)
        return detect_box_box(a, b, c);
    if (a->col_type==PHYS_COL_SPHERE && b->col_type==PHYS_COL_CAPSULE)
        return detect_sphere_capsule(a, b, c);
    if (a->col_type==PHYS_COL_CAPSULE&& b->col_type==PHYS_COL_SPHERE) {
        int r = detect_sphere_capsule(b, a, c);
        c->normal = pv3_neg(c->normal); return r;
    }
    if (a->col_type==PHYS_COL_CAPSULE&& b->col_type==PHYS_COL_CAPSULE)
        return detect_capsule_capsule(a, b, c);
    if (a->col_type==PHYS_COL_CAPSULE&& b->col_type==PHYS_COL_BOX)
        return detect_capsule_box(a, b, c);
    if (a->col_type==PHYS_COL_BOX    && b->col_type==PHYS_COL_CAPSULE) {
        int r = detect_capsule_box(b, a, c);
        c->normal = pv3_neg(c->normal); return r;
    }
    return 0;
}
/* ═══════════════════════════════════════════════════════
 *  phys.c  Part 3 — 衝突応答 / 動力学積分 / ジョイント
 * ═══════════════════════════════════════════════════════ */

/* ─── 衝突応答 (インパルスベース + 摩擦) ────────────────
 *  参考: Erin Catto GDC2006, Baumgarte stabilization
 * ──────────────────────────────────────────────────────── */
static void resolve_contact(PhysWorld* w, PhysContact* c) {
    PhysBody* a = &w->bodies[c->a];
    PhysBody* b = &w->bodies[c->b];

    /* トリガーは応答しない */
    if (a->is_trigger || b->is_trigger) return;
    /* 静的同士は不要 */
    if (a->body_type==PHYS_BODY_STATIC && b->body_type==PHYS_BODY_STATIC) return;

    float inv_ma = a->body_type==PHYS_BODY_DYNAMIC ? a->inv_mass : 0.0f;
    float inv_mb = b->body_type==PHYS_BODY_DYNAMIC ? b->inv_mass : 0.0f;

    float restitution = (a->restitution + b->restitution) * 0.5f;
    float friction    = (a->friction    + b->friction)    * 0.5f;

    /* Baumgarte 貫通補正 */
    float bias_factor  = 0.2f;
    float slop         = 0.005f;
    float bias = bias_factor * (c->depth - slop > 0 ? c->depth - slop : 0.0f);

    /* 各接触点でインパルス計算 */
    for (int i = 0; i < c->n_contacts; i++) {
        PhysVec3 n = c->normal;
        PhysVec3 ra = pv3_sub(c->contact[i], a->pos);
        PhysVec3 rb = pv3_sub(c->contact[i], b->pos);

        /* ワールド慣性テンソル (対角のみ近似) */
        PhysVec3 iia = a->body_type==PHYS_BODY_DYNAMIC ? a->inv_inertia : (PhysVec3){0,0,0};
        PhysVec3 iib = b->body_type==PHYS_BODY_DYNAMIC ? b->inv_inertia : (PhysVec3){0,0,0};

        /* 相対速度 */
        PhysVec3 va_at = pv3_add(a->vel, pv3_cross(a->omega, ra));
        PhysVec3 vb_at = pv3_add(b->vel, pv3_cross(b->omega, rb));
        PhysVec3 rel_vel = pv3_sub(va_at, vb_at);
        float vn = pv3_dot(rel_vel, n);

        /* 分離中なら何もしない */
        if (vn > 0.05f) continue;

        /* ─ 法線インパルス ─ */
        PhysVec3 ra_x_n = pv3_cross(ra, n);
        PhysVec3 rb_x_n = pv3_cross(rb, n);
        /* 有効質量の逆数 */
        float inv_eff = inv_ma + inv_mb
          + pv3_dot(ra_x_n, (PhysVec3){ra_x_n.x*iia.x, ra_x_n.y*iia.y, ra_x_n.z*iia.z})
          + pv3_dot(rb_x_n, (PhysVec3){rb_x_n.x*iib.x, rb_x_n.y*iib.y, rb_x_n.z*iib.z});
        if (inv_eff < P_EPS) continue;

        float jn = (-(1.0f + restitution) * vn + bias) / inv_eff;
        if (jn < 0.0f) jn = 0.0f;
        PhysVec3 impulse_n = pv3_scale(n, jn);

        /* 法線インパルス適用 */
        if (a->body_type==PHYS_BODY_DYNAMIC) {
            a->vel   = pv3_add(a->vel,   pv3_scale(impulse_n, inv_ma));
            PhysVec3 ang_a = pq_rot(a->rot, (PhysVec3){
                pv3_cross(ra,impulse_n).x * iia.x,
                pv3_cross(ra,impulse_n).y * iia.y,
                pv3_cross(ra,impulse_n).z * iia.z
            });
            a->omega = pv3_add(a->omega, ang_a);
        }
        if (b->body_type==PHYS_BODY_DYNAMIC) {
            b->vel   = pv3_sub(b->vel,   pv3_scale(impulse_n, inv_mb));
            PhysVec3 ang_b = pq_rot(b->rot, (PhysVec3){
                pv3_cross(rb,impulse_n).x * iib.x,
                pv3_cross(rb,impulse_n).y * iib.y,
                pv3_cross(rb,impulse_n).z * iib.z
            });
            b->omega = pv3_sub(b->omega, ang_b);
        }

        /* ─ 摩擦インパルス ─ */
        /* 再計算した相対速度 */
        va_at = pv3_add(a->vel, pv3_cross(a->omega, ra));
        vb_at = pv3_add(b->vel, pv3_cross(b->omega, rb));
        rel_vel = pv3_sub(va_at, vb_at);
        /* 接線方向 */
        PhysVec3 tangent = pv3_sub(rel_vel, pv3_scale(n, pv3_dot(rel_vel, n)));
        float tlen = pv3_len(tangent);
        if (tlen > P_EPS) {
            tangent = pv3_scale(tangent, 1.0f/tlen);
            float vt = pv3_dot(rel_vel, tangent);
            PhysVec3 ra_x_t = pv3_cross(ra, tangent);
            PhysVec3 rb_x_t = pv3_cross(rb, tangent);
            float inv_eff_t = inv_ma + inv_mb
              + pv3_dot(ra_x_t,(PhysVec3){ra_x_t.x*iia.x,ra_x_t.y*iia.y,ra_x_t.z*iia.z})
              + pv3_dot(rb_x_t,(PhysVec3){rb_x_t.x*iib.x,rb_x_t.y*iib.y,rb_x_t.z*iib.z});
            if (inv_eff_t > P_EPS) {
                float jt = -vt / inv_eff_t;
                /* コーロン摩擦 */
                float max_fric = friction * jn;
                jt = fclampf(jt, -max_fric, max_fric);
                PhysVec3 impulse_t = pv3_scale(tangent, jt);
                if (a->body_type==PHYS_BODY_DYNAMIC) {
                    a->vel   = pv3_add(a->vel,   pv3_scale(impulse_t,  inv_ma));
                    a->omega = pv3_add(a->omega,
                        pq_rot(a->rot,(PhysVec3){pv3_cross(ra,impulse_t).x*iia.x,
                                                  pv3_cross(ra,impulse_t).y*iia.y,
                                                  pv3_cross(ra,impulse_t).z*iia.z}));
                }
                if (b->body_type==PHYS_BODY_DYNAMIC) {
                    b->vel   = pv3_sub(b->vel,   pv3_scale(impulse_t,  inv_mb));
                    b->omega = pv3_sub(b->omega,
                        pq_rot(b->rot,(PhysVec3){pv3_cross(rb,impulse_t).x*iib.x,
                                                  pv3_cross(rb,impulse_t).y*iib.y,
                                                  pv3_cross(rb,impulse_t).z*iib.z}));
                }
            }
        }
    }
}

/* ─── 位置補正 (スペキュラー法) ─────────────────────── */
static void positional_correction(PhysWorld* w, PhysContact* c) {
    if (c->is_trigger) return;
    PhysBody* a = &w->bodies[c->a];
    PhysBody* b = &w->bodies[c->b];
    if (a->is_trigger || b->is_trigger) return;
    float inv_ma = a->body_type==PHYS_BODY_DYNAMIC ? a->inv_mass : 0.0f;
    float inv_mb = b->body_type==PHYS_BODY_DYNAMIC ? b->inv_mass : 0.0f;
    float total  = inv_ma + inv_mb;
    if (total < P_EPS) return;
    float correction = fmaxf(c->depth - 0.005f, 0.0f) / total * 0.4f;
    PhysVec3 cv = pv3_scale(c->normal, correction);
    if (a->body_type==PHYS_BODY_DYNAMIC) {
        if (!(a->freeze_pos&1)) a->pos.x -= cv.x * inv_ma;
        if (!(a->freeze_pos&2)) a->pos.y -= cv.y * inv_ma;
        if (!(a->freeze_pos&4)) a->pos.z -= cv.z * inv_ma;
    }
    if (b->body_type==PHYS_BODY_DYNAMIC) {
        if (!(b->freeze_pos&1)) b->pos.x += cv.x * inv_mb;
        if (!(b->freeze_pos&2)) b->pos.y += cv.y * inv_mb;
        if (!(b->freeze_pos&4)) b->pos.z += cv.z * inv_mb;
    }
}

/* ─── ジョイント解決 ─────────────────────────────────── */
static void solve_joints(PhysWorld* w, float dt) {
    for (int ji = 0; ji < PHYS_MAX_JOINTS; ji++) {
        PhysJoint* j = &w->joints[ji];
        if (!j->active) continue;
        if (j->a<0||j->a>=PHYS_MAX_BODIES||!w->bodies[j->a].active) continue;
        if (j->b<0||j->b>=PHYS_MAX_BODIES||!w->bodies[j->b].active) continue;
        PhysBody* ba = &w->bodies[j->a];
        PhysBody* bb = &w->bodies[j->b];
        float inv_ma = ba->body_type==PHYS_BODY_DYNAMIC ? ba->inv_mass : 0.0f;
        float inv_mb = bb->body_type==PHYS_BODY_DYNAMIC ? bb->inv_mass : 0.0f;

        switch (j->type) {
        case PHYS_JOINT_FIXED: {
            /* 位置制約: 相対位置を0に保つ */
            PhysVec3 ra = pq_rot(ba->rot, j->anchor_a);
            PhysVec3 rb = pq_rot(bb->rot, j->anchor_b);
            PhysVec3 pa = pv3_add(ba->pos, ra);
            PhysVec3 pb = pv3_add(bb->pos, rb);
            PhysVec3 err = pv3_sub(pb, pa);
            float total = inv_ma + inv_mb;
            if (total < P_EPS) break;
            float k = 0.5f / (dt * dt * total + 0.01f);
            PhysVec3 impulse = pv3_scale(err, k * dt);
            if (ba->body_type==PHYS_BODY_DYNAMIC)
                ba->vel = pv3_add(ba->vel, pv3_scale(impulse, inv_ma));
            if (bb->body_type==PHYS_BODY_DYNAMIC)
                bb->vel = pv3_sub(bb->vel, pv3_scale(impulse, inv_mb));
            break;
        }
        case PHYS_JOINT_DISTANCE: {
            PhysVec3 d = pv3_sub(bb->pos, ba->pos);
            float dist = pv3_len(d);
            if (dist < P_EPS) break;
            float error = 0.0f;
            if (dist < j->min_dist) error = dist - j->min_dist;
            else if (dist > j->max_dist) error = dist - j->max_dist;
            else break;
            PhysVec3 n = pv3_scale(d, 1.0f/dist);
            float total = inv_ma + inv_mb;
            if (total < P_EPS) break;
            float jn = error / total * 0.4f;
            PhysVec3 impulse = pv3_scale(n, jn);
            if (ba->body_type==PHYS_BODY_DYNAMIC)
                ba->pos = pv3_add(ba->pos, pv3_scale(impulse, inv_ma));
            if (bb->body_type==PHYS_BODY_DYNAMIC)
                bb->pos = pv3_sub(bb->pos, pv3_scale(impulse, inv_mb));
            break;
        }
        case PHYS_JOINT_SPRING: {
            PhysVec3 d = pv3_sub(bb->pos, ba->pos);
            float dist = pv3_len(d);
            if (dist < P_EPS) break;
            PhysVec3 n = pv3_scale(d, 1.0f/dist);
            /* バネ減衰力 */
            float force = -j->spring_k * (dist - j->rest_length);
            PhysVec3 rel_vel = pv3_sub(bb->vel, ba->vel);
            float damping = j->spring_damper * pv3_dot(rel_vel, n);
            float total_f = force - damping;
            PhysVec3 fvec = pv3_scale(n, total_f * dt);
            if (ba->body_type==PHYS_BODY_DYNAMIC)
                ba->vel = pv3_sub(ba->vel, pv3_scale(fvec, inv_ma));
            if (bb->body_type==PHYS_BODY_DYNAMIC)
                bb->vel = pv3_add(bb->vel, pv3_scale(fvec, inv_mb));
            break;
        }
        case PHYS_JOINT_HINGE: {
            /* ヒンジ: 角度制限のみ (簡易版) */
            PhysVec3 axis_w = pq_rot(ba->rot, j->axis);
            float rel_omega = pv3_dot(pv3_sub(bb->omega, ba->omega), axis_w);
            /* 角度は累積困難なので角速度のみ制限 */
            (void)rel_omega;
            /* 垂直方向の角速度を制約 */
            PhysVec3 perp = pv3_sub(bb->omega, pv3_scale(axis_w, pv3_dot(bb->omega, axis_w)));
            if (bb->body_type==PHYS_BODY_DYNAMIC)
                bb->omega = pv3_sub(bb->omega, pv3_scale(perp, 0.5f));
            break;
        }
        }
    }
}

/* ─── 積分 (Semi-implicit Euler) ─────────────────────── */
static void integrate_body(PhysBody* b, PhysVec3 gravity, float dt) {
    if (b->body_type != PHYS_BODY_DYNAMIC) return;
    if (b->is_sleeping) return;

    /* 重力 + 蓄積力を加速度に変換 */
    PhysVec3 accel = pv3_add(gravity,
                      pv3_scale(b->force, b->inv_mass));

    /* 線形速度 */
    b->vel = pv3_add(b->vel, pv3_scale(accel, dt));
    b->vel = pv3_scale(b->vel, fmaxf(0.0f, 1.0f - b->lin_drag * dt));

    /* 軸固定 */
    if (b->freeze_pos&1) b->vel.x = 0;
    if (b->freeze_pos&2) b->vel.y = 0;
    if (b->freeze_pos&4) b->vel.z = 0;

    /* 角速度 */
    PhysVec3 ang_accel = (PhysVec3){
        b->torque.x * b->inv_inertia.x,
        b->torque.y * b->inv_inertia.y,
        b->torque.z * b->inv_inertia.z
    };
    b->omega = pv3_add(b->omega, pv3_scale(ang_accel, dt));
    b->omega = pv3_scale(b->omega, fmaxf(0.0f, 1.0f - b->ang_drag * dt));
    if (b->freeze_rot&1) b->omega.x = 0;
    if (b->freeze_rot&2) b->omega.y = 0;
    if (b->freeze_rot&4) b->omega.z = 0;

    /* 位置更新 */
    if (!(b->freeze_pos&1)) b->pos.x += b->vel.x * dt;
    if (!(b->freeze_pos&2)) b->pos.y += b->vel.y * dt;
    if (!(b->freeze_pos&4)) b->pos.z += b->vel.z * dt;

    /* 回転更新 */
    b->rot = pq_integrate(b->rot, b->omega, dt);

    /* 力リセット */
    b->force  = (PhysVec3){0,0,0};
    b->torque = (PhysVec3){0,0,0};
}

/* ─── スリープ判定 ───────────────────────────────────── */
static void update_sleep(PhysBody* b, float thresh) {
    if (b->body_type != PHYS_BODY_DYNAMIC) return;
    float speed  = pv3_len(b->vel);
    float ospeed = pv3_len(b->omega);
    if (speed < thresh && ospeed < thresh) {
        b->sleep_counter++;
        if (b->sleep_counter >= PHYS_SLEEP_FRAMES) {
            b->is_sleeping   = 1;
            b->vel   = (PhysVec3){0,0,0};
            b->omega = (PhysVec3){0,0,0};
        }
    } else {
        b->sleep_counter = 0;
        b->is_sleeping   = 0;
    }
}

/* ─── レイヤーフィルタ ───────────────────────────────── */
static int layers_collide(PhysWorld* w, PhysBody* a, PhysBody* b) {
    if (!(a->layer_mask & (1u << b->layer))) return 0;
    if (!(b->layer_mask & (1u << a->layer))) return 0;
    if (!(w->layer_matrix[a->layer] & (1u << b->layer))) return 0;
    return 1;
}

/* ─── ワールドステップ ───────────────────────────────── */
void phys_world_step(PhysWorld* w, float dt) {
    if (!w || dt <= 0.0f) return;
    float sub_dt = dt / (float)w->substeps;

    for (int step = 0; step < w->substeps; step++) {
        /* 1. AABB 更新 */
        for (int i = 0; i < PHYS_MAX_BODIES; i++) {
            if (!w->bodies[i].active || w->bodies[i].col_type==PHYS_COL_NONE) continue;
            compute_aabb(&w->bodies[i]);
        }
        /* 2. 積分 */
        for (int i = 0; i < PHYS_MAX_BODIES; i++) {
            if (!w->bodies[i].active) continue;
            /* キネマティックは速度で移動 */
            if (w->bodies[i].body_type == PHYS_BODY_KINEMATIC) {
                w->bodies[i].pos = pv3_add(w->bodies[i].pos,
                                            pv3_scale(w->bodies[i].vel, sub_dt));
                w->bodies[i].rot = pq_integrate(w->bodies[i].rot,
                                                  w->bodies[i].omega, sub_dt);
                continue;
            }
            integrate_body(&w->bodies[i], w->gravity, sub_dt);
        }
        /* 3. 衝突検出 */
        w->n_contacts = 0;
        w->n_triggers = 0;
        for (int i = 0; i < PHYS_MAX_BODIES-1; i++) {
            if (!w->bodies[i].active||w->bodies[i].col_type==PHYS_COL_NONE) continue;
            for (int j = i+1; j < PHYS_MAX_BODIES; j++) {
                if (!w->bodies[j].active||w->bodies[j].col_type==PHYS_COL_NONE) continue;
                if (w->bodies[i].body_type==PHYS_BODY_STATIC &&
                    w->bodies[j].body_type==PHYS_BODY_STATIC) continue;
                if (!layers_collide(w, &w->bodies[i], &w->bodies[j])) continue;
                /* スリープ同士はスキップ */
                if (w->bodies[i].is_sleeping && w->bodies[j].is_sleeping) continue;
                /* ブロードフェーズ */
                if (!aabb_overlap(&w->bodies[i], &w->bodies[j])) continue;
                /* ナロウフェーズ */
                if (w->n_contacts >= PHYS_MAX_CONTACTS) break;
                PhysContact c;
                memset(&c, 0, sizeof(c));
                c.a = i; c.b = j;
                if (detect_pair(&w->bodies[i], &w->bodies[j], &c)) {
                    c.is_trigger = (w->bodies[i].is_trigger || w->bodies[j].is_trigger);
                    w->contacts[w->n_contacts++] = c;
                    if (c.is_trigger && w->n_triggers < PHYS_MAX_CONTACTS) {
                        w->trigger_pairs[w->n_triggers][0] = i;
                        w->trigger_pairs[w->n_triggers][1] = j;
                        w->n_triggers++;
                    }
                }
            }
        }
        /* 4. 衝突応答 */
        int iters = 6;
        for (int it = 0; it < iters; it++) {
            for (int ci = 0; ci < w->n_contacts; ci++)
                resolve_contact(w, &w->contacts[ci]);
        }
        /* 5. 位置補正 */
        for (int ci = 0; ci < w->n_contacts; ci++)
            positional_correction(w, &w->contacts[ci]);
        /* 6. ジョイント */
        solve_joints(w, sub_dt);
        /* 7. スリープ */
        for (int i = 0; i < PHYS_MAX_BODIES; i++) {
            if (w->bodies[i].active)
                update_sleep(&w->bodies[i], w->sleep_thresh);
        }
    }
}
/* ═══════════════════════════════════════════════════════
 *  phys.c  Part 4 — レイキャスト / オーバーラップ / クエリ
 * ═══════════════════════════════════════════════════════ */

/* ─── Ray vs Sphere ──────────────────────────────────── */
static int ray_sphere(PhysVec3 o, PhysVec3 d, PhysVec3 center, float r, float* t) {
    PhysVec3 oc = pv3_sub(o, center);
    float b = pv3_dot(oc, d);
    float c = pv3_dot(oc, oc) - r*r;
    float disc = b*b - c;
    if (disc < 0) return 0;
    float sq = sqrtf(disc);
    float t0 = -b - sq;
    float t1 = -b + sq;
    if (t1 < 0) return 0;
    *t = (t0 >= 0) ? t0 : t1;
    return 1;
}

/* ─── Ray vs AABB (slab) ─────────────────────────────── */
static int ray_aabb(PhysVec3 o, PhysVec3 inv_d, PhysVec3 mn, PhysVec3 mx, float max_t, float* t) {
    float tx1 = (mn.x - o.x)*inv_d.x, tx2 = (mx.x - o.x)*inv_d.x;
    float ty1 = (mn.y - o.y)*inv_d.y, ty2 = (mx.y - o.y)*inv_d.y;
    float tz1 = (mn.z - o.z)*inv_d.z, tz2 = (mx.z - o.z)*inv_d.z;
    float tmin = fmaxf(fmaxf(fminf(tx1,tx2), fminf(ty1,ty2)), fminf(tz1,tz2));
    float tmax = fminf(fminf(fmaxf(tx1,tx2), fmaxf(ty1,ty2)), fmaxf(tz1,tz2));
    if (tmax < 0 || tmin > tmax || tmin > max_t) return 0;
    *t = (tmin >= 0) ? tmin : tmax;
    return 1;
}

/* ─── Ray vs Capsule ─────────────────────────────────── */
static int ray_capsule(PhysVec3 ro, PhysVec3 rd,
                        PhysVec3 p0, PhysVec3 p1, float r, float* tout) {
    PhysVec3 ba = pv3_sub(p1, p0);
    PhysVec3 oa = pv3_sub(ro, p0);
    float baba = pv3_dot(ba,ba);
    float bard = pv3_dot(ba,rd);
    float baoa = pv3_dot(ba,oa);
    float rdoa = pv3_dot(rd,oa);
    float oaoa = pv3_dot(oa,oa);
    float a = baba - bard*bard;
    float b = baba*rdoa - baoa*bard;
    float c = baba*oaoa - baoa*baoa - r*r*baba;
    float h = b*b - a*c;
    if (h >= 0.0f) {
        float t2 = (-b - sqrtf(h))/a;
        float y2 = baoa + t2*bard;
        if (y2>0 && y2<baba) { if(t2>=0){ *tout=t2; return 1; } }
    }
    /* 端の球 */
    PhysVec3 oc = (h<0)?oa:pv3_sub(ro,p1);
    float boc = pv3_dot(ba,oc)/baba;
    if (boc<0) boc=0; else if (boc>1) boc=1;
    (void)boc;
    float t_sp;
    if (ray_sphere(ro,rd, p0, r, &t_sp) && t_sp>=0 && t_sp < *tout) { *tout=t_sp; return 1; }
    if (ray_sphere(ro,rd, p1, r, &t_sp) && t_sp>=0 && t_sp < *tout) { *tout=t_sp; return 1; }
    return 0;
}

/* ─── Ray vs OBB ─────────────────────────────────────── */
static int ray_obb(PhysVec3 ro, PhysVec3 rd,
                   PhysVec3 center, PhysQuat rot,
                   float hx, float hy, float hz, float max_t, float* t) {
    /* レイをOBBローカル空間へ変換 */
    PhysQuat inv = pq_conj(rot);
    PhysVec3 lo  = pq_rot(inv, pv3_sub(ro, center));
    PhysVec3 ld  = pq_rot(inv, rd);
    PhysVec3 inv_d = {
        fabsf(ld.x)>P_EPS ? 1.0f/ld.x : 1e30f,
        fabsf(ld.y)>P_EPS ? 1.0f/ld.y : 1e30f,
        fabsf(ld.z)>P_EPS ? 1.0f/ld.z : 1e30f
    };
    PhysVec3 mn = {-hx,-hy,-hz}, mx = {hx,hy,hz};
    return ray_aabb(lo, inv_d, mn, mx, max_t, t);
}

/* ─── 法線計算 ───────────────────────────────────────── */
static PhysVec3 hit_normal_sphere(PhysVec3 hit, PhysVec3 center) {
    return pv3_norm(pv3_sub(hit, center));
}
static PhysVec3 hit_normal_obb(PhysVec3 hit, PhysVec3 center, PhysQuat rot,
                                float hx, float hy, float hz) {
    PhysQuat inv = pq_conj(rot);
    PhysVec3 lp = pq_rot(inv, pv3_sub(hit, center));
    PhysVec3 ln = {0,1,0};
    float bx=fabsf(lp.x)-hx, by=fabsf(lp.y)-hy, bz=fabsf(lp.z)-hz;
    if (bx>by && bx>bz) ln=(PhysVec3){lp.x>0?1.0f:-1.0f,0,0};
    else if (by>bz)     ln=(PhysVec3){0,lp.y>0?1.0f:-1.0f,0};
    else                ln=(PhysVec3){0,0,lp.z>0?1.0f:-1.0f};
    return pq_rot(rot, ln);
}

/* ─── phys_raycast ───────────────────────────────────── */
PhysRayHit phys_raycast(PhysWorld* w, float ox, float oy, float oz,
                          float dx, float dy, float dz,
                          float max_dist, uint32_t layer_mask) {
    PhysRayHit ret = {0}; ret.t = max_dist;
    if (!w) return ret;
    PhysVec3 ro = {ox,oy,oz};
    PhysVec3 rd = pv3_norm((PhysVec3){dx,dy,dz});
    PhysVec3 inv_d = {
        fabsf(rd.x)>P_EPS?1.0f/rd.x:1e30f,
        fabsf(rd.y)>P_EPS?1.0f/rd.y:1e30f,
        fabsf(rd.z)>P_EPS?1.0f/rd.z:1e30f
    };

    for (int i = 0; i < PHYS_MAX_BODIES; i++) {
        PhysBody* b = &w->bodies[i];
        if (!b->active || b->col_type==PHYS_COL_NONE) continue;
        if (!(b->layer_mask & layer_mask)) continue;
        /* AABB 高速却下 */
        float ta;
        if (!ray_aabb(ro, inv_d, b->aabb_min, b->aabb_max, ret.t, &ta)) continue;
        /* 精密判定 */
        float t2;
        int hit = 0;
        PhysVec3 center = pv3_add(b->pos, pq_rot(b->rot, b->col_offset));
        switch (b->col_type) {
        case PHYS_COL_SPHERE:
            t2 = ret.t;
            hit = ray_sphere(ro, rd, center, b->col_size.x, &t2) && t2<ret.t && t2>=0;
            break;
        case PHYS_COL_BOX:
            t2 = ret.t;
            hit = ray_obb(ro, rd, center, b->rot, b->col_size.x, b->col_size.y, b->col_size.z, ret.t, &t2) && t2>=0;
            break;
        case PHYS_COL_CAPSULE: {
            PhysVec3 p0,p1; capsule_endpoints(b,&p0,&p1);
            t2 = ret.t;
            hit = ray_capsule(ro, rd, p0, p1, b->col_size.x, &t2) && t2>=0 && t2<ret.t;
            break;
        }
        default: break;
        }
        if (hit) {
            ret.hit     = 1;
            ret.body_id = i;
            ret.t       = t2;
            ret.point   = pv3_add(ro, pv3_scale(rd, t2));
            switch (b->col_type) {
            case PHYS_COL_SPHERE: ret.normal = hit_normal_sphere(ret.point, center); break;
            case PHYS_COL_BOX:    ret.normal = hit_normal_obb(ret.point, center, b->rot,
                                    b->col_size.x, b->col_size.y, b->col_size.z); break;
            case PHYS_COL_CAPSULE: ret.normal = hit_normal_sphere(ret.point,
                                    /* 最近傍端を使う近似 */
                                    pv3_lerp(pq_rot(b->rot,(PhysVec3){0,b->col_size.y*0.5f,0}),
                                             pq_rot(b->rot,(PhysVec3){0,-b->col_size.y*0.5f,0}),
                                             seg_closest_t(ret.point,
                                               pv3_add(b->pos,pq_rot(b->rot,(PhysVec3){0,b->col_size.y*0.5f,0})),
                                               pv3_add(b->pos,pq_rot(b->rot,(PhysVec3){0,-b->col_size.y*0.5f,0})))
                                            )); break;
            default: break;
            }
        }
    }
    return ret;
}

/* ─── phys_sphere_cast ───────────────────────────────── */
PhysRayHit phys_sphere_cast(PhysWorld* w,
                              float ox, float oy, float oz,
                              float dx, float dy, float dz,
                              float radius, float max_dist, uint32_t layer_mask) {
    /* 一時的な球ボディを作って各ボディと拡張済み判定 */
    PhysRayHit ret = {0}; ret.t = max_dist;
    if (!w) return ret;
    PhysVec3 ro = {ox,oy,oz};
    PhysVec3 rd = pv3_norm((PhysVec3){dx,dy,dz});

    for (int i = 0; i < PHYS_MAX_BODIES; i++) {
        PhysBody* b = &w->bodies[i];
        if (!b->active || b->col_type==PHYS_COL_NONE) continue;
        if (!(b->layer_mask & layer_mask)) continue;
        PhysVec3 center = pv3_add(b->pos, pq_rot(b->rot, b->col_offset));
        float expanded_r = b->col_size.x + radius;
        float t2 = ret.t;
        int hit = 0;
        /* 拡大した球でレイキャスト近似 */
        if (b->col_type == PHYS_COL_SPHERE) {
            hit = ray_sphere(ro, rd, center, expanded_r, &t2) && t2>=0 && t2<ret.t;
        } else {
            /* 拡大AABBで近似 */
            PhysVec3 mn = pv3_sub(b->aabb_min, (PhysVec3){radius,radius,radius});
            PhysVec3 mx = pv3_add(b->aabb_max, (PhysVec3){radius,radius,radius});
            PhysVec3 inv_d2 = {
                fabsf(rd.x)>P_EPS?1.0f/rd.x:1e30f,
                fabsf(rd.y)>P_EPS?1.0f/rd.y:1e30f,
                fabsf(rd.z)>P_EPS?1.0f/rd.z:1e30f
            };
            hit = ray_aabb(ro, inv_d2, mn, mx, ret.t, &t2) && t2>=0;
        }
        if (hit) {
            ret.hit = 1; ret.body_id = i; ret.t = t2;
            ret.point  = pv3_add(ro, pv3_scale(rd, t2));
            ret.normal = hit_normal_sphere(ret.point, center);
        }
    }
    return ret;
}

/* ─── phys_overlap_sphere ────────────────────────────── */
int phys_overlap_sphere(PhysWorld* w, float cx, float cy, float cz,
                         float radius, uint32_t layer_mask, int* out, int max_out) {
    if (!w || !out) return 0;
    int count = 0;
    PhysVec3 c = {cx,cy,cz};
    for (int i = 0; i < PHYS_MAX_BODIES && count < max_out; i++) {
        PhysBody* b = &w->bodies[i];
        if (!b->active || b->col_type==PHYS_COL_NONE) continue;
        if (!(b->layer_mask & layer_mask)) continue;
        /* 中心距離で簡易判定 */
        PhysVec3 bc = pv3_add(b->pos, pq_rot(b->rot, b->col_offset));
        float dist2 = pv3_len2(pv3_sub(c, bc));
        float r_sum = radius + b->col_size.x;
        if (dist2 <= r_sum*r_sum) out[count++] = i;
    }
    return count;
}

/* ─── phys_overlap_box ───────────────────────────────── */
int phys_overlap_box(PhysWorld* w, float cx, float cy, float cz,
                      float hx, float hy, float hz,
                      uint32_t layer_mask, int* out, int max_out) {
    if (!w || !out) return 0;
    int count = 0;
    PhysVec3 qmn={(cx-hx),(cy-hy),(cz-hz)}, qmx={(cx+hx),(cy+hy),(cz+hz)};
    for (int i = 0; i < PHYS_MAX_BODIES && count < max_out; i++) {
        PhysBody* b = &w->bodies[i];
        if (!b->active || b->col_type==PHYS_COL_NONE) continue;
        if (!(b->layer_mask & layer_mask)) continue;
        if (b->aabb_min.x<=qmx.x && b->aabb_max.x>=qmn.x &&
            b->aabb_min.y<=qmx.y && b->aabb_max.y>=qmn.y &&
            b->aabb_min.z<=qmx.z && b->aabb_max.z>=qmn.z)
            out[count++] = i;
    }
    return count;
}

/* ─── 衝突/トリガー確認 ──────────────────────────────── */
int phys_check_collision(PhysWorld* w, int a, int b) {
    if (!w) return 0;
    for (int i = 0; i < w->n_contacts; i++) {
        if (!w->contacts[i].is_trigger &&
            ((w->contacts[i].a==a && w->contacts[i].b==b) ||
             (w->contacts[i].a==b && w->contacts[i].b==a)))
            return 1;
    }
    return 0;
}
int phys_check_trigger(PhysWorld* w, int a, int b) {
    if (!w) return 0;
    for (int i = 0; i < w->n_triggers; i++) {
        if ((w->trigger_pairs[i][0]==a && w->trigger_pairs[i][1]==b) ||
            (w->trigger_pairs[i][0]==b && w->trigger_pairs[i][1]==a))
            return 1;
    }
    return 0;
}
