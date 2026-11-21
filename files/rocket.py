import math, numpy as np, cadquery as cq

# ---------- parameters (mm) ----------
OD, WALL = 26.0, 1.0            # airframe
BODY_L = 200.0
MOTOR_D, MOTOR_L = 18.0, 70.0   # standard 18 mm black-powder motor
SLEEVE_ID, SLEEVE_WALL = 18.4, 0.8
BLOCK_Z = 64.0                  # thrust ring height -> motor sticks out 6 mm
NOSE_L, SHOULDER_L, FIT = 70.0, 15.0, 0.1
N_FINS, ROOT, TIP, SPAN, SWEEP, FIN_T = 3, 60.0, 30.0, 48.0, 30.0, 2.0
LUG_ID, LUG_WALL, LUG_L, LUG_Z = 4.4, 1.0, 30.0, 60.0
R, RI = OD / 2, OD / 2 - WALL

# ---------- airframe: tube + motor mount + fins + launch lug ----------
tube = cq.Workplane("XY").circle(R).circle(RI).extrude(BODY_L)
rs_i, rs_o = SLEEVE_ID / 2, SLEEVE_ID / 2 + SLEEVE_WALL
sleeve = cq.Workplane("XY").circle(rs_o).circle(rs_i).extrude(BLOCK_Z)
ring_aft = cq.Workplane("XY").circle(RI + 0.1).circle(rs_o - 0.1).extrude(2)
ring_fwd = ring_aft.translate((0, 0, BLOCK_Z - 2))
# thrust ring with 45-degree underside so it prints without support
block = (cq.Workplane("XZ")
         .polyline([(rs_i, BLOCK_Z), (RI + 0.1, BLOCK_Z), (RI + 0.1, BLOCK_Z + 5),
                    (6.5, BLOCK_Z + 5), (6.5, BLOCK_Z + rs_i - 6.5)]).close()
         .revolve(360, (0, 0, 0), (0, 1, 0)))
fin = (cq.Workplane("XZ")
       .polyline([(R - 0.5, 0), (R - 0.5, ROOT), (R + SPAN, ROOT - SWEEP), (R + SPAN, ROOT - SWEEP - TIP)]).close()
       .extrude(FIN_T / 2, both=True))
airframe = tube.union(sleeve).union(ring_aft).union(ring_fwd).union(block)
for k in range(N_FINS):
    airframe = airframe.union(fin.rotate((0, 0, 0), (0, 0, 1), k * 360 / N_FINS))
lug_r = R + LUG_ID / 2 + 0.4
lug = (cq.Workplane("XY").circle(LUG_ID / 2 + LUG_WALL).extrude(LUG_L)
       .translate((lug_r, 0, LUG_Z)).rotate((0, 0, 0), (0, 0, 1), 60))
lug_hole = (cq.Workplane("XY").circle(LUG_ID / 2).extrude(LUG_L + 2)
            .translate((lug_r, 0, LUG_Z - 1)).rotate((0, 0, 0), (0, 0, 1), 60))
inner = cq.Workplane("XY").circle(RI).extrude(BODY_L)
airframe = airframe.union(lug.cut(inner)).cut(lug_hole)


# ---------- nose cone: hollow tangent ogive + shoulder + shock-cord bar ----------
def ogive(r, l):
    rho = (r * r + l * l) / (2 * r)
    zm = l / 2
    xm = math.sqrt(rho * rho - zm * zm) + r - rho
    return (cq.Workplane("XZ").moveTo(0, 0).lineTo(r, 0).threePointArc((xm, zm), (0, l)).close()
            .revolve(360, (0, 0, 0), (0, 1, 0)))


rsh = RI - FIT
nose = ogive(R, NOSE_L).cut(ogive(R - 1.2, NOSE_L - 6))
shoulder = cq.Workplane("XY").workplane(offset=-SHOULDER_L).circle(rsh).circle(rsh - 1.2).extrude(SHOULDER_L)
bar = (cq.Workplane("XY").box(2 * rsh, 3, 3, centered=(True, True, False)).translate((0, 0, -SHOULDER_L))
       .intersect(cq.Workplane("XY").workplane(offset=-SHOULDER_L).circle(rsh).extrude(3)))
nose = nose.union(shoulder).union(bar)

# ---------- export ----------
for name, part in [("rocket_airframe", airframe), ("rocket_nosecone", nose)]:
    cq.exporters.export(part, f"{name}.stl", tolerance=0.02, angularTolerance=0.1)
    cq.exporters.export(part, f"{name}.step")
assy = cq.Assembly()
assy.add(airframe, name="airframe", color=cq.Color(0.9, 0.9, 0.92))
assy.add(nose, name="nosecone", loc=cq.Location(cq.Vector(0, 0, BODY_L)), color=cq.Color(0.85, 0.15, 0.1))
assy.save("rocket_assembly.step")

# ---------- mass, CG, Barrowman CP ----------
RHO_PLA = 1.24e-3                       # g/mm3, assumes walls print solid
MOTOR_MASS = 24.0                       # g, approx loaded C6-5 (assumption)
TOTAL = NOSE_L + BODY_L
va, vn = airframe.val().Volume(), nose.val().Volume()
ca = cq.Shape.centerOfMass(airframe.val()).z
cn = cq.Shape.centerOfMass(nose.val()).z + BODY_L
ma, mn = va * RHO_PLA, vn * RHO_PLA
motor_z = MOTOR_L / 2 - (MOTOR_L - BLOCK_Z)
NOSE_WEIGHT = 10.0                      # g of clay pressed into the nose tip
nw_z = BODY_L + NOSE_L - 30
m_tot = ma + mn + MOTOR_MASS + NOSE_WEIGHT
cg_z = (ma * ca + mn * cn + MOTOR_MASS * motor_z + NOSE_WEIGHT * nw_z) / m_tot
cg = TOTAL - cg_z                       # station from nose tip

cn_nose, x_nose = 2.0, 0.466 * NOSE_L
lf = math.hypot(SPAN, SWEEP + TIP / 2 - ROOT / 2)
cn_fin = (1 + R / (SPAN + R)) * (4 * N_FINS * (SPAN / OD) ** 2) / (1 + math.sqrt(1 + (2 * lf / (ROOT + TIP)) ** 2))
x_fin = (TOTAL - ROOT) + SWEEP / 3 * (ROOT + 2 * TIP) / (ROOT + TIP) + (ROOT + TIP - ROOT * TIP / (ROOT + TIP)) / 6
cp = (cn_nose * x_nose + cn_fin * x_fin) / (cn_nose + cn_fin)
print(f"airframe {ma:.1f} g | nose {mn:.1f} g | with motor + nose weight {m_tot:.1f} g")
print(f"length {TOTAL:.0f} mm | CG {cg:.1f} mm | CP {cp:.1f} mm from nose tip")
print(f"static margin {(cp - cg) / OD:.2f} calibers")

# ---------- preview: outside view + cutaway ----------
import trimesh, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

half = cq.Workplane("XY").box(200, 200, 400, centered=(False, True, True))   # keep x >= 0
views = []
for cut in (False, True):
    pair = []
    for part, dz in ((airframe, 0), (nose, BODY_L)):
        p = part.intersect(half) if cut else part
        cq.exporters.export(p, "_tmp.stl", tolerance=0.05, angularTolerance=0.2)
        m = trimesh.load("_tmp.stl"); m.apply_translation([0, 0, dz]); pair.append(m)
    views.append(pair)
light = np.array([0.3, -0.6, 0.7]); light /= np.linalg.norm(light)
from PIL import Image, ImageChops
tiles = []
for pair, title in zip(views, ["assembled", "cutaway: motor mount, thrust ring, hollow nose"]):
    fig = plt.figure(figsize=(16, 16)); ax = fig.add_subplot(projection="3d")
    for mesh, c in zip(pair, ([0.85, 0.87, 0.9], [0.85, 0.2, 0.12])):
        sh = np.abs(mesh.face_normals[:, [2, 0, 1]] @ light) * 0.7 + 0.3
        ax.add_collection3d(Poly3DCollection(mesh.triangles[:, :, [2, 0, 1]], facecolors=np.array(c) * sh[:, None], edgecolors="none"))
    ax.set_xlim(0, TOTAL); ax.set_ylim(-62, 62); ax.set_zlim(-62, 62)
    ax.set_box_aspect((TOTAL, 124, 124)); ax.view_init(12, -80); ax.set_axis_off()
    fig.savefig("_v.png", dpi=110); plt.close(fig)
    im = Image.open("_v.png").convert("RGB")
    im = im.crop(ImageChops.difference(im, Image.new("RGB", im.size, "white")).getbbox())
    fig = plt.figure(figsize=(im.width / 100, im.height / 100 + 0.5)); plt.imshow(im); plt.axis("off"); plt.title(title, fontsize=13)
    fig.savefig("_v.png", dpi=100, bbox_inches="tight"); plt.close(fig); tiles.append(Image.open("_v.png").convert("RGB"))
W = max(t.width for t in tiles); out = Image.new("RGB", (W, sum(t.height for t in tiles) + 20), "white"); y = 0
for t in tiles:
    out.paste(t, ((W - t.width) // 2, y)); y += t.height + 20
out.save("rocket_preview.png")
