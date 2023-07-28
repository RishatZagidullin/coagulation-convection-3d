#version 330 core
layout (location = 0) in vec3 aPos;

uniform mat4 mvp;
uniform vec4 clip;

void main()
{
    gl_Position = mvp*vec4(aPos, 1.0);
    gl_ClipDistance[0] = dot(clip, vec4(aPos, 1.0));
}
