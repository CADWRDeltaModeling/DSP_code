stream_points = [(540750.0394097868, 4184077.0644629044), (545820.8043179004, 4185853.8432327043), (551293.782989001, 4188427.8341326755), (552076.5000000003, 4190882.8000000026), (549540.8999999976, 4193652.800000003), (548812.5999999997, 4196773.400000003), (548940.0999999994, 4201215.500000001), (550154.0999999996, 4204018.7), (554728.699999999, 4209338.800000003), (563174.2999999995, 4212397.500000002), (565861.1000000001, 4213416.299999999), (571461.6000000002, 4212373.0), (572901.999999999, 4210177.699999997), (574896.5999999989, 4209435.000000002), (581863.7570720065, 4213066.606876341), (584078.9000000003, 4213518.399999999), (586759.7350863256, 4213024.954723995), (590874.9000000026, 4212727.199999999), (595023.4000000014, 4211296.700000001), (597949.7000000014, 4211961.8999999985), (600340.2999999999, 4213411.2), (603382.4000000019, 4213340.8), (606090.600000001, 4213609.0), (615408.5, 4224528.6000000015)]

elev_clip = 1

for i in range ( len ( stream_points ) -1):
    AddPlot("Pseudocolor", "salinity")
    p = PseudocolorAttributes()
    p.minFlag = 1
    p.min = 0.0
    p.maxFlag = 1
    p.max = 35.0
    p.colorTableName = "hot"
    p.legendFlag =0
    SetPlotOptions(p)
    AddOperator("Transform")
    tr = TransformAttributes()
    tr.doScale = 1
    tr.scaleZ = 600
    SetOperatorOptions(tr)
    x1 = stream_points[i][0]
    y1 = stream_points[i][1]
    x2 = stream_points[i+1][0]
    y2 = stream_points[i+1][1]
    normal0 = y1 - y2
    normal1 = x2 - x1
    AddOperator("Slice")
    s = SliceAttributes()
    s.originType = s.Point
    s.project2d = 0
    s.originPoint =(x1, y1, 0)
    s.normal = (normal0, normal1, 0)
    SetOperatorOptions(s)
    AddOperator("Clip")
    c = ClipAttributes()
    c.plane1Origin = (x1, y1, 0)
    c.plane1Normal = (1, 0, 0)
    c.plane2Origin =(x2, y2, 0)
    c.plane2Normal = (-1, 0, 0)
    c.planeInverse = 1
    c.plane2Status = 1
    SetOperatorOptions(c)
    DrawPlots()


v = View3DAttributes()
v.viewNormal = (0.114452, -0.650205, 0.751088)
v.focus = (540730, 4.22455e+06, 1338.23)
v.viewUp = (0.0903594, 0.759736, 0.643923)
v.viewAngle = 30
v.parallelScale = 109730
v.nearPlane = -219460
v.farPlane = 219460
v.imagePan = (-0.141668, 0.00512045)
v.imageZoom = 2.09868
v.perspective = 1
v.eyeAngle = 2
v.shear = (0, 0, 1)
SetView3D(v) # Set the 3D view
# AddPlot("Pseudocolor", "element_level")
# p = PseudocolorAttributes()
# p.legendFlag = 0
# SetPlotOptions(p)
# AddOperator("Elevate")
# DrawPlots()
