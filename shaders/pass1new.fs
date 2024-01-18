#version 330 core

in vec3 pix;
out vec4 fragColor;

// ----------------------------------------------------------------------------- //

uniform uint frameCounter;
uniform int maxBounce;
uniform int nTriangles;
uniform int nNodes;
uniform int width;
uniform int height;
uniform int hdrResolution;

uniform samplerBuffer triangles;
uniform samplerBuffer nodes;

uniform sampler2D lastFrame;
uniform sampler2D hdrMap;
uniform sampler2D hdrCache;

uniform vec3 eye;
uniform mat4 cameraRotate;

// ----------------------------------------------------------------------------- //

#define PI              3.1415926
#define INF             114514.0
#define SIZE_TRIANGLE   12
#define SIZE_BVHNODE    4

// ----------------------------------------------------------------------------- //

// Triangle ���ݸ�ʽ
struct Triangle {
    vec3 p1, p2, p3;    // ��������
    vec3 n1, n2, n3;    // ���㷨��
};

// BVH ���ڵ�
struct BVHNode {
    int left;           // ������
    int right;          // ������
    int n;              // ������������Ŀ
    int index;          // ����������
    vec3 AA, BB;        // ��ײ��
};

// ���������ʶ���
struct Material {
    vec3 emissive;          // ��Ϊ��Դʱ�ķ�����ɫ
    vec3 baseColor;
    float subsurface;
    float metallic;
    float specular;
    float specularTint;
    float roughness;
    float anisotropic;
    float sheen;
    float sheenTint;
    float clearcoat;
    float clearcoatGloss;
    float IOR;
    float transmission;
};

// ����
struct Ray {
    vec3 startPoint;
    vec3 direction;
};

// �����󽻽��
struct HitResult {
    bool isHit;             // �Ƿ�����
    bool isInside;          // �Ƿ���ڲ�����
    float distance;         // �뽻��ľ���
    vec3 hitPoint;          // �������е�
    vec3 normal;            // ���е㷨��
    vec3 viewDir;           // ���иõ�Ĺ��ߵķ���
    Material material;      // ���е�ı������
};

// ��Ҫ�Բ����ķ��ؽ��
struct SampleResult {
    vec3 direction;
    float pdf;
};

// ----------------------------------------------------------------------------- //

// ��ȡ�� i �±��������
Triangle getTriangle(int i) {
    int offset = i * SIZE_TRIANGLE;
    Triangle t;

    // ��������
    t.p1 = texelFetch(triangles, offset + 0).xyz;
    t.p2 = texelFetch(triangles, offset + 1).xyz;
    t.p3 = texelFetch(triangles, offset + 2).xyz;
    // ����
    t.n1 = texelFetch(triangles, offset + 3).xyz;
    t.n2 = texelFetch(triangles, offset + 4).xyz;
    t.n3 = texelFetch(triangles, offset + 5).xyz;

    return t;
}

// ��ȡ�� i �±�������εĲ���
Material getMaterial(int i) {
    Material m;

    int offset = i * SIZE_TRIANGLE;
    vec3 param1 = texelFetch(triangles, offset + 8).xyz;
    vec3 param2 = texelFetch(triangles, offset + 9).xyz;
    vec3 param3 = texelFetch(triangles, offset + 10).xyz;
    vec3 param4 = texelFetch(triangles, offset + 11).xyz;
    
    m.emissive = texelFetch(triangles, offset + 6).xyz;
    m.baseColor = texelFetch(triangles, offset + 7).xyz;
    m.subsurface = param1.x;
    m.metallic = param1.y;
    m.specular = param1.z;
    m.specularTint = param2.x;
    m.roughness = param2.y;
    m.anisotropic = param2.z;
    m.sheen = param3.x;
    m.sheenTint = param3.y;
    m.clearcoat = param3.z;
    m.clearcoatGloss = param4.x;
    m.IOR = param4.y;
    m.transmission = param4.z;

    return m;
}

// ��ȡ�� i �±�� BVHNode ����
BVHNode getBVHNode(int i) {
    BVHNode node;

    // ��������
    int offset = i * SIZE_BVHNODE;
    ivec3 childs = ivec3(texelFetch(nodes, offset + 0).xyz);
    ivec3 leafInfo = ivec3(texelFetch(nodes, offset + 1).xyz);
    node.left = int(childs.x);
    node.right = int(childs.y);
    node.n = int(leafInfo.x);
    node.index = int(leafInfo.y);

    // ��Χ��
    node.AA = texelFetch(nodes, offset + 2).xyz;
    node.BB = texelFetch(nodes, offset + 3).xyz;

    return node;
}

// ----------------------------------------------------------------------------- //

// ���ߺ��������� 
HitResult hitTriangle(Triangle triangle, Ray ray) {
    HitResult res;
    res.distance = INF;
    res.isHit = false;
    res.isInside = false;

    vec3 p1 = triangle.p1;
    vec3 p2 = triangle.p2;
    vec3 p3 = triangle.p3;

    vec3 S = ray.startPoint;    // �������
    vec3 d = ray.direction;     // ���߷���
    vec3 N = normalize(cross(p2-p1, p3-p1));    // ������

    // �������α���ģ���ڲ�������
    if (dot(N, d) > 0.0f) {
        N = -N;   
        res.isInside = true;
    }

    // ������ߺ�������ƽ��
    if (abs(dot(N, d)) < 0.00001f) return res;

    // ����
    float t = (dot(N, p1) - dot(S, N)) / dot(d, N);
    if (t < 0.0005f) return res;    // ����������ڹ��߱���

    // �������
    vec3 P = S + d * t;

    // �жϽ����Ƿ�����������
    vec3 c1 = cross(p2 - p1, P - p1);
    vec3 c2 = cross(p3 - p2, P - p2);
    vec3 c3 = cross(p1 - p3, P - p3);
    bool r1 = (dot(c1, N) > 0 && dot(c2, N) > 0 && dot(c3, N) > 0);
    bool r2 = (dot(c1, N) < 0 && dot(c2, N) < 0 && dot(c3, N) < 0);

    // ���У���װ���ؽ��
    if (r1 || r2) {
        res.isHit = true;
        res.hitPoint = P;
        res.distance = t;
        res.normal = N;
        res.viewDir = d;
        
        // ���ݽ���λ�ò�ֵ���㷨��
        float alpha = (-(P.x-p2.x)*(p3.y-p2.y) + (P.y-p2.y)*(p3.x-p2.x)) / (-(p1.x-p2.x)*(p3.y-p2.y) + (p1.y-p2.y)*(p3.x-p2.x)+1e-7);
        float beta  = (-(P.x-p3.x)*(p1.y-p3.y) + (P.y-p3.y)*(p1.x-p3.x)) / (-(p2.x-p3.x)*(p1.y-p3.y) + (p2.y-p3.y)*(p1.x-p3.x)+1e-7);
        float gama  = 1.0 - alpha - beta;
        vec3 Nsmooth = alpha * triangle.n1 + beta * triangle.n2 + gama * triangle.n3;
        Nsmooth = normalize(Nsmooth);
        /*
        vec3 Nsmooth = N;*/
        res.normal = (res.isInside) ? (-Nsmooth) : (Nsmooth);
    }

    return res;
}

// �� aabb �����󽻣�û�н����򷵻� -1
float hitAABB(Ray r, vec3 AA, vec3 BB) {
    vec3 invdir = 1.0 / r.direction;

    vec3 f = (BB - r.startPoint) * invdir;
    vec3 n = (AA - r.startPoint) * invdir;

    vec3 tmax = max(f, n);
    vec3 tmin = min(f, n);

    float t1 = min(tmax.x, min(tmax.y, tmax.z));
    float t0 = max(tmin.x, max(tmin.y, tmin.z));

    return (t1 >= t0) ? ((t0 > 0.0) ? (t0) : (t1)) : (-1);
}

// ----------------------------------------------------------------------------- //

// �������������±귶Χ [l, r] ���������
HitResult hitArray(Ray ray, int l, int r) {
    HitResult res;
    res.isHit = false;
    res.distance = INF;
    for(int i=l; i<=r; i++) {
        Triangle triangle = getTriangle(i);
        HitResult r = hitTriangle(triangle, ray);
        if(r.isHit && r.distance<res.distance) {
            res = r;
            res.material = getMaterial(i);
        }
    }
    return res;
}

// ���� BVH ��
HitResult hitBVH(Ray ray) {
    HitResult res;
    res.isHit = false;
    res.distance = INF;

    // ջ
    int stack[256];
    int sp = 0;

    stack[sp++] = 1;
    while(sp>0) {
        int top = stack[--sp];
        BVHNode node = getBVHNode(top);
        
        // ��Ҷ�ӽڵ㣬���������Σ����������
        if(node.n>0) {
            int L = node.index;
            int R = node.index + node.n - 1;
            HitResult r = hitArray(ray, L, R);
            if(r.isHit && r.distance<res.distance) res = r;
            continue;
        }
        
        // �����Һ��� AABB ��
        float d1 = INF; // ����Ӿ���
        float d2 = INF; // �Һ��Ӿ���
        if(node.left>0) {
            BVHNode leftNode = getBVHNode(node.left);
            d1 = hitAABB(ray, leftNode.AA, leftNode.BB);
        }
        if(node.right>0) {
            BVHNode rightNode = getBVHNode(node.right);
            d2 = hitAABB(ray, rightNode.AA, rightNode.BB);
        }

        // ������ĺ���������
        if(d1>0 && d2>0) {
            if(d1<d2) { // d1<d2, �����
                stack[sp++] = node.right;
                stack[sp++] = node.left;
            } else {    // d2<d1, �ұ���
                stack[sp++] = node.left;
                stack[sp++] = node.right;
            }
        } else if(d1>0) {   // ���������
            stack[sp++] = node.left;
        } else if(d2>0) {   // �������ұ�
            stack[sp++] = node.right;
        }
    }

    return res;
}

// ----------------------------------------------------------------------------- //


uint seed = uint(
    uint((pix.x * 0.5 + 0.5) * width)  * uint(1973) + 
    uint((pix.y * 0.5 + 0.5) * height) * uint(9277) + 
    uint(frameCounter) * uint(26699)) | uint(1);

uint wang_hash(inout uint seed) {
    seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
    seed *= uint(9);
    seed = seed ^ (seed >> 4);
    seed *= uint(0x27d4eb2d);
    seed = seed ^ (seed >> 15);
    return seed;
}
 
float rand() {
    return float(wang_hash(seed)) / 4294967296.0;
}

uint seed_sync = uint(
    uint((pix.x * 0.0 + 0.5) * width)  * uint(1973) + 
    uint((pix.y * 0.0 + 0.5) * height) * uint(9277) + 
    uint(114514) * uint(26699)) | uint(1);

float rand_sync() {
    return float(wang_hash(seed_sync)) / 4294967296.0;
}

// ----------------------------------------------------------------------------- //

// 1 ~ 8 ά�� sobol ���ɾ���
const uint V[8*32] = uint[8*32](
    uint(2147483648), uint(1073741824), uint(536870912), uint(268435456), uint(134217728), uint(67108864), uint(33554432), uint(16777216), uint(8388608), uint(4194304), uint(2097152), uint(1048576), uint(524288), uint(262144), uint(131072), uint(65536), uint(32768), uint(16384), uint(8192), uint(4096), uint(2048), uint(1024), uint(512), uint(256), uint(128), uint(64), uint(32), uint(16), uint(8), uint(4), uint(2), uint(1),
    uint(2147483648), uint(3221225472), uint(2684354560), uint(4026531840), uint(2281701376), uint(3422552064), uint(2852126720), uint(4278190080), uint(2155872256), uint(3233808384), uint(2694840320), uint(4042260480), uint(2290614272), uint(3435921408), uint(2863267840), uint(4294901760), uint(2147516416), uint(3221274624), uint(2684395520), uint(4026593280), uint(2281736192), uint(3422604288), uint(2852170240), uint(4278255360), uint(2155905152), uint(3233857728), uint(2694881440), uint(4042322160), uint(2290649224), uint(3435973836), uint(2863311530), uint(4294967295),
    uint(2147483648), uint(3221225472), uint(1610612736), uint(2415919104), uint(3892314112), uint(1543503872), uint(2382364672), uint(3305111552), uint(1753219072), uint(2629828608), uint(3999268864), uint(1435500544), uint(2154299392), uint(3231449088), uint(1626210304), uint(2421489664), uint(3900735488), uint(1556135936), uint(2388680704), uint(3314585600), uint(1751705600), uint(2627492864), uint(4008611328), uint(1431684352), uint(2147543168), uint(3221249216), uint(1610649184), uint(2415969680), uint(3892340840), uint(1543543964), uint(2382425838), uint(3305133397),
    uint(2147483648), uint(3221225472), uint(536870912), uint(1342177280), uint(4160749568), uint(1946157056), uint(2717908992), uint(2466250752), uint(3632267264), uint(624951296), uint(1507852288), uint(3872391168), uint(2013790208), uint(3020685312), uint(2181169152), uint(3271884800), uint(546275328), uint(1363623936), uint(4226424832), uint(1977167872), uint(2693105664), uint(2437829632), uint(3689389568), uint(635137280), uint(1484783744), uint(3846176960), uint(2044723232), uint(3067084880), uint(2148008184), uint(3222012020), uint(537002146), uint(1342505107),
    uint(2147483648), uint(1073741824), uint(536870912), uint(2952790016), uint(4160749568), uint(3690987520), uint(2046820352), uint(2634022912), uint(1518338048), uint(801112064), uint(2707423232), uint(4038066176), uint(3666345984), uint(1875116032), uint(2170683392), uint(1085997056), uint(579305472), uint(3016343552), uint(4217741312), uint(3719483392), uint(2013407232), uint(2617981952), uint(1510979072), uint(755882752), uint(2726789248), uint(4090085440), uint(3680870432), uint(1840435376), uint(2147625208), uint(1074478300), uint(537900666), uint(2953698205),
    uint(2147483648), uint(1073741824), uint(1610612736), uint(805306368), uint(2818572288), uint(335544320), uint(2113929216), uint(3472883712), uint(2290089984), uint(3829399552), uint(3059744768), uint(1127219200), uint(3089629184), uint(4199809024), uint(3567124480), uint(1891565568), uint(394297344), uint(3988799488), uint(920674304), uint(4193267712), uint(2950604800), uint(3977188352), uint(3250028032), uint(129093376), uint(2231568512), uint(2963678272), uint(4281226848), uint(432124720), uint(803643432), uint(1633613396), uint(2672665246), uint(3170194367),
    uint(2147483648), uint(3221225472), uint(2684354560), uint(3489660928), uint(1476395008), uint(2483027968), uint(1040187392), uint(3808428032), uint(3196059648), uint(599785472), uint(505413632), uint(4077912064), uint(1182269440), uint(1736704000), uint(2017853440), uint(2221342720), uint(3329785856), uint(2810494976), uint(3628507136), uint(1416089600), uint(2658719744), uint(864310272), uint(3863387648), uint(3076993792), uint(553150080), uint(272922560), uint(4167467040), uint(1148698640), uint(1719673080), uint(2009075780), uint(2149644390), uint(3222291575),
    uint(2147483648), uint(1073741824), uint(2684354560), uint(1342177280), uint(2281701376), uint(1946157056), uint(436207616), uint(2566914048), uint(2625634304), uint(3208642560), uint(2720006144), uint(2098200576), uint(111673344), uint(2354315264), uint(3464626176), uint(4027383808), uint(2886631424), uint(3770826752), uint(1691164672), uint(3357462528), uint(1993345024), uint(3752330240), uint(873073152), uint(2870150400), uint(1700563072), uint(87021376), uint(1097028000), uint(1222351248), uint(1560027592), uint(2977959924), uint(23268898), uint(437609937)
);

// ������ 
uint grayCode(uint i) {
	return i ^ (i>>1);
}

// ���ɵ� d ά�ȵĵ� i �� sobol ��
float sobol(uint d, uint i) {
    uint result = uint(0);
    uint offset = d * uint(32);
    for(uint j = uint(0); bool(i); i >>= uint(1), j++) 
        if(bool(i & uint(1)))
            result ^= V[j+offset];

    return float(result) * (1.0f/float(0xFFFFFFFFU));
}

// ���ɵ� i ֡�ĵ� b �η�����Ҫ�Ķ�ά�������
vec2 sobolVec2(uint i, uint b) {
    float u = sobol(b*uint(2), grayCode(i));
    float v = sobol(b*uint(2)+uint(1), grayCode(i));
    return vec2(u, v);
}

vec2 CranleyPattersonRotation(vec2 p) {
    uint pseed = uint(
        uint((pix.x * 0.5 + 0.5) * width)  * uint(1973) + 
        uint((pix.y * 0.5 + 0.5) * height) * uint(9277) + 
        uint(114514/1919) * uint(26699)) | uint(1);
    
    float u = float(wang_hash(pseed)) / 4294967296.0;
    float v = float(wang_hash(pseed)) / 4294967296.0;

    p.x += u;
    if(p.x>1) p.x -= 1;
    if(p.x<0) p.x += 1;

    p.y += v;
    if(p.y>1) p.y -= 1;
    if(p.y<0) p.y += 1;

    return p;
}

// ----------------------------------------------------------------------------- //

float sqr(float x) { 
    return x*x; 
}

float SchlickFresnel(float u) {
    float m = clamp(1-u, 0, 1);
    float m2 = m*m;
    return m2*m2*m; // pow(m,5)
}

float GTR1(float NdotH, float a) {
    if (a >= 1) return 1/PI;
    float a2 = a*a;
    float t = 1 + (a2-1)*NdotH*NdotH;
    return (a2-1) / (PI*log(a2)*t);
}

float GTR2(float NdotH, float a) {
    float a2 = a*a;
    float t = 1 + (a2-1)*NdotH*NdotH;
    return a2 / (PI * t*t);
}

float GTR2_aniso(float NdotH, float HdotX, float HdotY, float ax, float ay) {
    return 1 / (PI * ax*ay * sqr( sqr(HdotX/ax) + sqr(HdotY/ay) + NdotH*NdotH ));
}

float smithG_GGX(float NdotV, float alphaG) {
    float a = alphaG*alphaG;
    float b = NdotV*NdotV;
    return 1 / (NdotV + sqrt(a + b - a*b));
}

float smithG_GGX_aniso(float NdotV, float VdotX, float VdotY, float ax, float ay) {
    return 1 / (NdotV + sqrt( sqr(VdotX*ax) + sqr(VdotY*ay) + sqr(NdotV) ));
}

vec3 BRDF_Evaluate_aniso(vec3 V, vec3 N, vec3 L, vec3 X, vec3 Y, in Material material) {
    float NdotL = dot(N, L);
    float NdotV = dot(N, V);
    if(NdotL < 0 || NdotV < 0) return vec3(0);

    vec3 H = normalize(L + V);
    float NdotH = dot(N, H);
    float LdotH = dot(L, H);

    // ������ɫ
    vec3 Cdlin = material.baseColor;
    float Cdlum = 0.3 * Cdlin.r + 0.6 * Cdlin.g  + 0.1 * Cdlin.b;
    vec3 Ctint = (Cdlum > 0) ? (Cdlin/Cdlum) : (vec3(1));   
    vec3 Cspec = material.specular * mix(vec3(1), Ctint, material.specularTint);
    vec3 Cspec0 = mix(0.08*Cspec, Cdlin, material.metallic); // 0�� ���淴����ɫ
    vec3 Csheen = mix(vec3(1), Ctint, material.sheenTint);   // ֯����ɫ

    // ������
    float Fd90 = 0.5 + 2.0 * LdotH * LdotH * material.roughness;
    float FL = SchlickFresnel(NdotL);
    float FV = SchlickFresnel(NdotV);
    float Fd = mix(1.0, Fd90, FL) * mix(1.0, Fd90, FV);

    // �α���ɢ��
    float Fss90 = LdotH * LdotH * material.roughness;
    float Fss = mix(1.0, Fss90, FL) * mix(1.0, Fss90, FV);
    float ss = 1.25 * (Fss * (1.0 / (NdotL + NdotV) - 0.5) + 0.5);
     
    // ���淴�� -- ����ͬ��
    float alpha = max(0.001, sqr(material.roughness));
    float Ds = GTR2(NdotH, alpha);
    float FH = SchlickFresnel(LdotH);
    vec3 Fs = mix(Cspec0, vec3(1), FH);
    float Gs = smithG_GGX(NdotL, material.roughness);
    Gs *= smithG_GGX(NdotV, material.roughness);
    /*
    // ���淴�� -- ��������
    float aspect = sqrt(1.0 - material.anisotropic * 0.9);
    float ax = max(0.001, sqr(material.roughness)/aspect);
    float ay = max(0.001, sqr(material.roughness)*aspect);
    float Ds = GTR2_aniso(NdotH, dot(H, X), dot(H, Y), ax, ay);
    float FH = SchlickFresnel(LdotH);
    vec3 Fs = mix(Cspec0, vec3(1), FH);
    float Gs;
    Gs  = smithG_GGX_aniso(NdotL, dot(L, X), dot(L, Y), ax, ay);
    Gs *= smithG_GGX_aniso(NdotV, dot(V, X), dot(V, Y), ax, ay);
    */

    // ����
    float Dr = GTR1(NdotH, mix(0.1, 0.001, material.clearcoatGloss));
    float Fr = mix(0.04, 1.0, FH);
    float Gr = smithG_GGX(NdotL, 0.25) * smithG_GGX(NdotV, 0.25);

    // sheen
    vec3 Fsheen = FH * material.sheen * Csheen;
    
    vec3 diffuse = (1.0/PI) * mix(Fd, ss, material.subsurface) * Cdlin + Fsheen;
    vec3 specular = Gs * Fs * Ds;
    vec3 clearcoat = vec3(0.25 * Gr * Fr * Dr * material.clearcoat);

    return diffuse * (1.0 - material.metallic) + specular + clearcoat;
}

vec3 BRDF_Evaluate(vec3 V, vec3 N, vec3 L, in Material material) {
    float NdotL = dot(N, L);
    float NdotV = dot(N, V);
    if(NdotL < 0 || NdotV < 0) return vec3(0);

    vec3 H = normalize(L + V);
    float NdotH = dot(N, H);
    float LdotH = dot(L, H);

    // ������ɫ
    vec3 Cdlin = material.baseColor;
    float Cdlum = 0.3 * Cdlin.r + 0.6 * Cdlin.g  + 0.1 * Cdlin.b;
    vec3 Ctint = (Cdlum > 0) ? (Cdlin/Cdlum) : (vec3(1));   
    vec3 Cspec = material.specular * mix(vec3(1), Ctint, material.specularTint);
    vec3 Cspec0 = mix(0.08*Cspec, Cdlin, material.metallic); // 0�� ���淴����ɫ
    vec3 Csheen = mix(vec3(1), Ctint, material.sheenTint);   // ֯����ɫ

    // ������
    float Fd90 = 0.5 + 2.0 * LdotH * LdotH * material.roughness;
    float FL = SchlickFresnel(NdotL);
    float FV = SchlickFresnel(NdotV);
    float Fd = mix(1.0, Fd90, FL) * mix(1.0, Fd90, FV);

    // �α���ɢ��
    float Fss90 = LdotH * LdotH * material.roughness;
    float Fss = mix(1.0, Fss90, FL) * mix(1.0, Fss90, FV);
    float ss = 1.25 * (Fss * (1.0 / (NdotL + NdotV) - 0.5) + 0.5);
     
    // ���淴�� -- ����ͬ��
    float alpha = max(0.001, sqr(material.roughness));
    float Ds = GTR2(NdotH, alpha);
    float FH = SchlickFresnel(LdotH);
    vec3 Fs = mix(Cspec0, vec3(1), FH);
    float Gs = smithG_GGX(NdotL, material.roughness);
    Gs *= smithG_GGX(NdotV, material.roughness);

    // ����
    float Dr = GTR1(NdotH, mix(0.1, 0.001, material.clearcoatGloss));
    float Fr = mix(0.04, 1.0, FH);
    float Gr = smithG_GGX(NdotL, 0.25) * smithG_GGX(NdotV, 0.25);

    // sheen
    vec3 Fsheen = FH * material.sheen * Csheen;
    
    vec3 diffuse = (1.0/PI) * mix(Fd, ss, material.subsurface) * Cdlin + Fsheen;
    vec3 specular = Gs * Fs * Ds;
    vec3 clearcoat = vec3(0.25 * Gr * Fr * Dr * material.clearcoat);

    return diffuse * (1.0 - material.metallic) + specular + clearcoat;
}

// ----------------------------------------------------------------------------- //

void getTangent(vec3 N, inout vec3 tangent, inout vec3 bitangent) {
    vec3 helper = vec3(1, 0, 0);
    if(abs(N.x)>0.999) helper = vec3(0, 0, 1);
    bitangent = normalize(cross(N, helper));
    tangent = normalize(cross(N, bitangent));
}

// ������ v ͶӰ�� N �ķ������
vec3 toNormalHemisphere(vec3 v, vec3 N) {
    vec3 helper = vec3(1, 0, 0);
    if(abs(N.x)>0.999) helper = vec3(0, 0, 1);
    vec3 tangent = normalize(cross(N, helper));
    vec3 bitangent = normalize(cross(N, tangent));
    return v.x * tangent + v.y * bitangent + v.z * N;
}

// ������Ȳ���
vec3 SampleHemisphere(float xi_1, float xi_2) {
    //xi_1 = rand(), xi_2 = rand();
    float z = xi_1;
    float r = max(0, sqrt(1.0 - z*z));
    float phi = 2.0 * PI * xi_2;
    return vec3(r * cos(phi), r * sin(phi), z);
}

// ���Ҽ�Ȩ�ķ���������
vec3 SampleCosineHemisphere(float xi_1, float xi_2, vec3 N) {
    // ���Ȳ��� xy Բ��Ȼ��ͶӰ�� z ����
    float r = sqrt(xi_1);
    float theta = xi_2 * 2.0 * PI;
    float x = r * cos(theta);
    float y = r * sin(theta);
    float z = sqrt(1.0 - x*x - y*y);

    // �� z ����ͶӰ���������
    vec3 L = toNormalHemisphere(vec3(x, y, z), N);
    return L;
}

// GTR2 ��Ҫ�Բ���
vec3 SampleGTR2(float xi_1, float xi_2, vec3 V, vec3 N, float alpha) {
    
    float phi_h = 2.0 * PI * xi_1;
    float sin_phi_h = sin(phi_h);
    float cos_phi_h = cos(phi_h);

    float cos_theta_h = sqrt((1.0-xi_2)/(1.0+(alpha*alpha-1.0)*xi_2));
    float sin_theta_h = sqrt(max(0.0, 1.0 - cos_theta_h * cos_theta_h));

    // ���� "΢ƽ��" �ķ����� ��Ϊ���淴��İ������ h 
    vec3 H = vec3(sin_theta_h*cos_phi_h, sin_theta_h*sin_phi_h, cos_theta_h);
    H = toNormalHemisphere(H, N);   // ͶӰ�������ķ������

    // ���� "΢����" ���㷴��ⷽ��
    vec3 L = reflect(-V, H);

    return L;
}

// GTR1 ��Ҫ�Բ���
vec3 SampleGTR1(float xi_1, float xi_2, vec3 V, vec3 N, float alpha) {
    
    float phi_h = 2.0 * PI * xi_1;
    float sin_phi_h = sin(phi_h);
    float cos_phi_h = cos(phi_h);

    float cos_theta_h = sqrt((1.0-pow(alpha*alpha, 1.0-xi_2))/(1.0-alpha*alpha));
    float sin_theta_h = sqrt(max(0.0, 1.0 - cos_theta_h * cos_theta_h));

    // ���� "΢ƽ��" �ķ����� ��Ϊ���淴��İ������ h 
    vec3 H = vec3(sin_theta_h*cos_phi_h, sin_theta_h*sin_phi_h, cos_theta_h);
    H = toNormalHemisphere(H, N);   // ͶӰ�������ķ������

    // ���� "΢����" ���㷴��ⷽ��
    vec3 L = reflect(-V, H);

    return L;
}

// ���շ���ȷֲ��ֱ�������� BRDF
vec3 SampleBRDF(float xi_1, float xi_2, float xi_3, vec3 V, vec3 N, in Material material) {
    float alpha_GTR1 = mix(0.1, 0.001, material.clearcoatGloss);
    float alpha_GTR2 = max(0.001, sqr(material.roughness));
    
    // �����ͳ��
    float r_diffuse = (1.0 - material.metallic);
    float r_specular = 1.0;
    float r_clearcoat = 0.25 * material.clearcoat;
    float r_sum = r_diffuse + r_specular + r_clearcoat;

    // ���ݷ���ȼ������
    float p_diffuse = r_diffuse / r_sum;
    float p_specular = r_specular / r_sum;
    float p_clearcoat = r_clearcoat / r_sum;

    // ���ո��ʲ���
    float rd = xi_3;

    // ������
    if(rd <= p_diffuse) {
        return SampleCosineHemisphere(xi_1, xi_2, N);
    } 
    // ���淴��
    else if(p_diffuse < rd && rd <= p_diffuse + p_specular) {    
        return SampleGTR2(xi_1, xi_2, V, N, alpha_GTR2);
    } 
    // ����
    else if(p_diffuse + p_specular < rd) {
        return SampleGTR1(xi_1, xi_2, V, N, alpha_GTR1);
    }
    return vec3(0, 1, 0);
}

// ����Ԥ����� HDR cache
vec3 SampleHdr(float xi_1, float xi_2) {
    vec2 xy = texture2D(hdrCache, vec2(xi_1, xi_2)).rg; // x, y
    xy.y = 1.0 - xy.y; // flip y

    // ��ȡ�Ƕ�
    float phi = 2.0 * PI * (xy.x - 0.5);    // [-pi ~ pi]
    float theta = PI * (xy.y - 0.5);        // [-pi/2 ~ pi/2]   

    // ��������㷽��
    vec3 L = vec3(cos(theta)*cos(phi), sin(theta), cos(theta)*sin(phi));

    return L;
}

// ----------------------------------------------------------------------------- //

// ����ά���� v תΪ HDR map ���������� uv
vec2 toSphericalCoord(vec3 v) {
    vec2 uv = vec2(atan(v.z, v.x), asin(v.y));
    uv /= vec2(2.0 * PI, PI);
    uv += 0.5;
    uv.y = 1.0 - uv.y;
    return uv;
}

// ��ȡ HDR ������ɫ
vec3 hdrColor(vec3 L) {
    vec2 uv = toSphericalCoord(normalize(L));
    vec3 color = texture2D(hdrMap, uv).rgb;
    return color;
}

// ������߷��� L ��ȡ HDR �ڸ�λ�õĸ����ܶ�
// hdr �ֱ���Ϊ 4096 x 2048 --> hdrResolution = 4096
float hdrPdf(vec3 L, int hdrResolution) {
    vec2 uv = toSphericalCoord(normalize(L));   // ��������ת uv ��������

    float pdf = texture2D(hdrCache, uv).b;      // ���������ܶ�
    float theta = PI * (0.5 - uv.y);            // theta ��Χ [-pi/2 ~ pi/2]
    float sin_theta = max(sin(theta), 1e-10);

    // �������ͼƬ�������ת��ϵ��
    float p_convert = float(hdrResolution * hdrResolution / 2) / (2.0 * PI * PI * sin_theta);  
    
    return pdf * p_convert;
}

// ��ȡ BRDF �� L �����ϵĸ����ܶ�
float BRDF_Pdf(vec3 V, vec3 N, vec3 L, in Material material) {
    float NdotL = dot(N, L);
    float NdotV = dot(N, V);
    if(NdotL < 0 || NdotV < 0) return 0;

    vec3 H = normalize(L + V);
    float NdotH = dot(N, H);
    float LdotH = dot(L, H);
     
    // ���淴�� -- ����ͬ��
    float alpha = max(0.001, sqr(material.roughness));
    float Ds = GTR2(NdotH, alpha); 
    float Dr = GTR1(NdotH, mix(0.1, 0.001, material.clearcoatGloss));   // ����

    // �ֱ�������� BRDF �ĸ����ܶ�
    float pdf_diffuse = NdotL / PI;
    float pdf_specular = Ds * NdotH / (4.0 * dot(L, H));
    float pdf_clearcoat = Dr * NdotH / (4.0 * dot(L, H));

    // �����ͳ��
    float r_diffuse = (1.0 - material.metallic);
    float r_specular = 1.0;
    float r_clearcoat = 0.25 * material.clearcoat;
    float r_sum = r_diffuse + r_specular + r_clearcoat;

    // ���ݷ���ȼ���ѡ��ĳ�ֲ�����ʽ�ĸ���
    float p_diffuse = r_diffuse / r_sum;
    float p_specular = r_specular / r_sum;
    float p_clearcoat = r_clearcoat / r_sum;

    // ���ݸ��ʻ�� pdf
    float pdf = p_diffuse   * pdf_diffuse 
              + p_specular  * pdf_specular
              + p_clearcoat * pdf_clearcoat;

    pdf = max(1e-10, pdf);
    return pdf;
}

float misMixWeight(float a, float b) {
    float t = a * a;
    return t / (b*b + t);
}

// ----------------------------------------------------------------------------- //

// ·��׷��
vec3 pathTracing(HitResult hit, int maxBounce) {

    vec3 Lo = vec3(0);      // ���յ���ɫ
    vec3 history = vec3(1); // �ݹ���۵���ɫ

    for(int bounce=0; bounce<maxBounce; bounce++) {
        vec3 V = -hit.viewDir;
        vec3 N = hit.normal;
        
        vec2 uv = sobolVec2(frameCounter+uint(1), uint(bounce));
        uv = CranleyPattersonRotation(uv);
        //uv = vec2(rand(), rand());

        vec3 L = SampleHemisphere(uv.x, uv.y);
        L = toNormalHemisphere(L, hit.normal);                          // ���䷽�� wi
        float pdf = 1.0 / (2.0 * PI);                                   // ������Ȳ��������ܶ�
        float cosine_o = max(0, dot(V, N));                             // �����ͷ��߼н�����
        float cosine_i = max(0, dot(L, N));                             // �����ͷ��߼н�����
        vec3 tangent, bitangent;
        getTangent(N, tangent, bitangent);
        vec3 f_r = BRDF_Evaluate_aniso(V, N, L, tangent, bitangent, hit.material);

        // �������
        Ray randomRay;
        randomRay.startPoint = hit.hitPoint;
        randomRay.direction = L;
        HitResult newHit = hitBVH(randomRay);

        // δ����
        if(!newHit.isHit) {
            vec3 skyColor = hdrColor(randomRay.direction);
            Lo += history * skyColor * f_r * cosine_i / pdf;
            break;
        }
        
        // ���й�Դ������ɫ
        vec3 Le = newHit.material.emissive;
        Lo += history * Le * f_r * cosine_i / pdf;
        
        // �ݹ�(����)
        hit = newHit;
        history *= f_r * cosine_i / pdf;  // �ۻ���ɫ
    }
    
    return Lo;
}

// ·��׷�� -- ��Ҫ�Բ����汾
vec3 pathTracingImportanceSampling(HitResult hit, int maxBounce) {

    vec3 Lo = vec3(0);      // ���յ���ɫ
    vec3 history = vec3(1); // �ݹ���۵���ɫ

    for(int bounce=0; bounce<maxBounce; bounce++) {
        vec3 V = -hit.viewDir;
        vec3 N = hit.normal;       

        // HDR ������ͼ��Ҫ�Բ���    
        Ray hdrTestRay;
        hdrTestRay.startPoint = hit.hitPoint;
        hdrTestRay.direction = SampleHdr(rand(), rand());

        // ����һ���󽻲��� �ж��Ƿ����ڵ�
        if(dot(N, hdrTestRay.direction) > 0.0) { // �������������� p ���������, ��Ϊ N dot L < 0            
            HitResult hdrHit = hitBVH(hdrTestRay);
            
            // ��չ����û���ڵ�������»�������
            if(!hdrHit.isHit) {
                // ��ȡ�������� L �ϵ�: 1.���չ���, 2.������ͼ�ڸ�λ�õ� pdf, 3.BRDF ����ֵ, 4.BRDF �ڸ÷���� pdf
                vec3 L = hdrTestRay.direction;
                vec3 color = hdrColor(L);
                float pdf_light = hdrPdf(L, hdrResolution);
                vec3 f_r = BRDF_Evaluate(V, N, L, hit.material);
                float pdf_brdf = BRDF_Pdf(V, N, L, hit.material);
                
                // ������Ҫ�Բ���
                float mis_weight = misMixWeight(pdf_light, pdf_brdf);
                Lo += mis_weight * history * color * f_r * dot(N, L) / pdf_light;
                //Lo += history * color * f_r * dot(N, L) / pdf_light;   // ����
            }
        }
        
        // ��ȡ 3 �������
        vec2 uv = sobolVec2(frameCounter+uint(1), uint(bounce));
        uv = CranleyPattersonRotation(uv);
        float xi_1 = uv.x;
        float xi_2 = uv.y;
        float xi_3 = rand();    // xi_3 �Ǿ��������������, ���� rand �ͺ�

        // ���� BRDF �õ�һ������ L
        vec3 L = SampleBRDF(xi_1, xi_2, xi_3, V, N, hit.material); 
        float NdotL = dot(N, L);
        if(NdotL <= 0.0) break;

        // �������
        Ray randomRay;
        randomRay.startPoint = hit.hitPoint;
        randomRay.direction = L;
        HitResult newHit = hitBVH(randomRay);

        // ��ȡ L �����ϵ� BRDF ֵ�͸����ܶ�
        vec3 f_r = BRDF_Evaluate(V, N, L, hit.material);
        float pdf_brdf = BRDF_Pdf(V, N, L, hit.material);
        if(pdf_brdf <= 0.0) break;

        // δ����        
        if(!newHit.isHit) {
            vec3 color = hdrColor(L);
            float pdf_light = hdrPdf(L, hdrResolution);            
            
            // ������Ҫ�Բ���
            float mis_weight = misMixWeight(pdf_brdf, pdf_light);   // f(a,b) = a^2 / (a^2 + b^2)
            Lo += mis_weight * history * color * f_r * NdotL / pdf_brdf;
            //Lo += history * color * f_r * NdotL / pdf_brdf;   // ����

            break;
        }
        
        // ���й�Դ������ɫ
        vec3 Le = newHit.material.emissive;
        Lo += history * Le * f_r * NdotL / pdf_brdf;             

        // �ݹ�(����)
        hit = newHit;
        history *= f_r * NdotL / pdf_brdf;   // �ۻ���ɫ
    }
    
    return Lo;
}

// ----------------------------------------------------------------------------- //

void main() {
    /*
    if(pix.y>0) {
        fragColor = texture2D(hdrMap, vec2(pix.x*0.5+0.5, 1-pix.y));
    } else {
        fragColor = texture2D(hdrCache, vec2(pix.x*0.5+0.5, -pix.y));
    }
    */

    // �������
    /*
    float pdf = texture2D(hdrCache, vec2(pix.x*0.5+0.5, -pix.y*0.5+0.5)).b;
    fragColor = vec4(pdf) * 100000;
    
    fragColor = texture2D(hdrMap, vec2(pix.x*0.5+0.5, -pix.y*0.5+0.5));

    // ������
    for(int i=0; i<500; i++) {
        vec2 uv =sobolVec2(i+1, 0); // ȡ 1,2 ά�ȵ� sobol ���������
        vec2 texcoord = texture2D(hdrCache, uv).rg;
        texcoord.y = 1.0 - texcoord.y;
        if(distance(pix.xy*0.5+0.5, texcoord)<0.005) fragColor.rgb = vec3(1, 0, 0);
    }*/
    
    /**/
    // Ͷ�����
    Ray ray;
    
    ray.startPoint = eye;
    vec2 AA = vec2((rand()-0.5)/float(width), (rand()-0.5)/float(height));
    vec4 dir = cameraRotate * vec4(pix.xy+AA, -1.5, 0.0);
    ray.direction = normalize(dir.xyz);

    // primary hit
    HitResult firstHit = hitBVH(ray);
    vec3 color;
    
    if(!firstHit.isHit) {
        color = vec3(0);
        color = hdrColor(ray.direction);
    } else {
        vec3 Le = firstHit.material.emissive;
        vec3 Li = pathTracingImportanceSampling(firstHit, maxBounce);
        //vec3 Li = pathTracing(firstHit, maxBounce);
        color = Le + Li;
    }
    
    // ����һ֡���
    vec3 lastColor = texture2D(lastFrame, pix.xy*0.5+0.5).rgb;
    color = mix(lastColor, color, 1.0/float(frameCounter+uint(1)));

    // ���
    gl_FragData[0] = vec4(color, 1.0);
    
}
