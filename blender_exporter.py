bl_info = {
    "name": "wlme",
    "description": "Writes geometry and specific attributes to disk that are used in the Wavelike Mapping Experiment.",
    "author": "Jeff Redbeard (jeffredbeard7@gmail.com)",
    "version": (1, 0),
    "blender": (2, 93, 0),
    "location": "File > Export > WLME",
    "warning": "",
    "wiki_url": "",
    "tracker_url": "",
    "support": 'COMMUNITY',
    "category": "Import-Export"
}

import bpy
import bmesh
import struct #https://docs.python.org/3/library/struct.html
from bpy import context
from bpy_extras.io_utils import ExportHelper
import json
   
def writeObject(self, context):
    ob = context.object
    uv_layer = ob.data.uv_layers.active.data
    
    vertBuff = []
    uvBuff   = []
    faceBuff = []
    #rebuild vertex, uv and face indices excluding duplicates
    for poly in ob.data.polygons:
        for index in poly.loop_indices:
            thisVertex = ob.data.vertices[ob.data.loops[index].vertex_index].co 
            thisUV = uv_layer[index].uv
            
            #check if already in the list
            i = 0
            found = 0
            for v in vertBuff:
                if(abs(v.x-thisVertex.x) <= max(1e-09 * max(abs(v.x), abs(thisVertex.x)), 0.0)):
                    if(abs(v.y-thisVertex.y) <= max(1e-09 * max(abs(v.y), abs(thisVertex.y)), 0.0)):
                        if(abs(v.z-thisVertex.z) <= max(1e-09 * max(abs(v.z), abs(thisVertex.z)), 0.0)):
                            if(abs(uvBuff[i].x-thisUV.x) <= max(1e-09 * max(abs(uvBuff[i].x), abs(thisUV.x)), 0.0)):
                                if(abs(uvBuff[i].y-thisUV.y) <= max(1e-09 * max(abs(uvBuff[i].y), abs(thisUV.y)), 0.0)):
                                    faceBuff.append(int(i))
                                    found = 1
                                    break
                i+=1
            #otherwise stash a new vertex
            if(found==0):
                faceBuff.append(len(vertBuff)) #index
                vertBuff.append(thisVertex)    #float, float, float
                uvBuff.append(thisUV)          #float, float
                
    #write to file
    if(self.format == "OPT_A"):
        with open(self.filepath, 'w') as ofile:
            ofile.write("%d " % len(vertBuff)) #num unique vertex/uv pairs
            ofile.write("%d " % len(faceBuff)) #num indices
            for v in vertBuff:
                ofile.write("%f %f %f " % v[:])
            for t in uvBuff:
                ofile.write("%f %f " % t[:])
            for p in faceBuff:
                ofile.write("%d " % p)
            ofile.close()
        return {'FINISHED'}
    else:
        
        with open(self.filepath, 'wb') as ofile:
            ofile.write(struct.pack('H', len(vertBuff))) 
            ofile.write(struct.pack('H', len(faceBuff)))
        
            for v in vertBuff:
                ofile.write(struct.pack('3f', v.x, v.y, v.z)) #v[:])) #"%f %f %f " % v[:])
            for t in uvBuff:
                ofile.write(struct.pack('2f', t.x, t.y)) #t[:])) #"%f %f " % t[:])
            for p in faceBuff:
                ofile.write(struct.pack('H', p)) #"%d " % p)
            ofile.close()
        return {'FINISHED'}    

class ObjectExport(bpy.types.Operator, ExportHelper):
    """My object export script"""
    bl_idname = "object.export_wlme"
    bl_label = "WLME Format Export"
    filename_ext = ".wlmadata"
    
    filter_glob     = bpy.props.StringProperty(default="*.wlmadata", options={'HIDDEN'}, maxlen=255)
    use_triangulate = bpy.props.BoolProperty(name="Triangulate", description="Triangulate object", default=True)

    def write_object_mesh(self, obj, fdj):
        out = []

        version = 1
        # This should be analogous to specularity.
        reflect = obj.get('reflect', 0.5)
        # This should be analogous to smoothness. With zero
        # there is no random variation in the reflection
        # vector. This unit is expected to be in degrees. 
        reflect_vector_noise = obj.get('reflect_vector_noise', 0.0) 

        print('obj', obj)
        print('obj.data', obj.data)
        
        world_matrix = obj.matrix_world
        data = bmesh.new()
        data.from_mesh(obj.data)
        data.transform(world_matrix)

        for face in data.faces:
            out_face = []
            for vert in face.verts:
                vertex = vert.co
                out_face.append([vertex[0], vertex[1], vertex[2]])
                print('writing vert')
            print('writing face')
            out.append(out_face)
        '''
          We can do a lot here. Using nodes we can produce a image surface
          that represents reflectivity and absorption. Then, finally, we can
          write that data out here, in binary if needed, which can be used
          by the simulator... That is a lot of work but it does become a big
          possibility with Blender doing the heavy lifting.
        '''
        fdj['meshes'].append({
          'name': obj.name,
          'polygons': out,
          'meta': [version, reflect, reflect_vector_noise],
        })

        data.free()

    def write_object_empty(self, obj, fdj):
        if 'wlme_type' in obj:
            wlme_type = obj.get('wlme_type', 'unknown')
        else:
            return
        if wlme_type == 'tx':
            wlme_power = obj.get('wlme_power', 1.0)
            fdj[wlme_type].append({
              'name': obj.name,
              'location': [obj.location[0], obj.location[1], obj.location[2]],
              'power': wlme_power,
            })
        elif wlme_type == 'rx':
            wlme_sensitivity = obj.get('wlme_sensitivity', 1.0)
            fdj[wlme_type].append({
              'name': obj.name,
              'location': [obj.location[0], obj.location[1], obj.location[2]],
              'sensitivity': wlme_sensitivity,
            })
        else:
            # Do not know what WLME type this is.
            print('unknown empty wlme_type field')
            return

    def triangulate_object(self, obj):
        print('triangulating object', obj)
        objdata = obj.data
        bm = bmesh.new()
        bm.from_mesh(objdata)
        bmesh.ops.triangulate(bm, faces=bm.faces[:]) #, quad_method=0, ngon_method=0)
        bm.to_mesh(objdata)
        bm.free()
         
    def execute(self, context):
        fdj = {
          'tx': [],
          'rx': [],
          'meshes': [],
        }
        # This will iterate through all the objects in the scene.
        objs = context.scene.objects
        for obj in objs:
            if obj.mode == 'EDIT':
                obj.mode_set(mode='OBJECT')
    
            print(obj.type)
            
            if obj.type == 'MESH':
                print('got mesh')
                # convert to triangles write to output
                self.triangulate_object(obj)
                self.write_object_mesh(obj, fdj)
            elif obj.type == 'SURFACE':
                # pull out nurbs convert to mesh then write to output
                self.triangulate_object(obj)
                self.write_object_mesh(obj, fdj)
            elif obj.type == 'EMPTY':
                print('got empty')
                self.write_object_empty(obj, fdj)
            fp = self.filepath
            print('fp=%s' % fp)
            with open(fp, 'w') as fd:
                fd.write(json.dumps(fdj))
        return {'FINISHED'}

    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}



# Add trigger into a dynamic menu
def menu_func_export(self, context):
    self.layout.operator(ObjectExport.bl_idname, text="WLME Export (.wlmedata)")
    

def register():
    bpy.utils.register_class(ObjectExport)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)


def unregister():
    bpy.utils.unregister_class(ObjectExport)
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)


if __name__ == "__main__":
    register()
