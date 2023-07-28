#pragma once
#include <GLFW/glfw3.h>
#include <glm/gtc/type_ptr.hpp>
#include "shader.h"

namespace gl
{
    enum drawable{grid, object, vector_field};

    class grid_viewer
    {
        private:
        glm::vec3 camera_pos;
        glm::vec4 clip;
        glm::mat4 mvp;

        Shader *fill_shader;
        Shader *line_shader;
        Shader *vect_shader;
        Shader *obj_shader;

        float aaa[16] = {
            1.810660, 0.000000, 0.000000, 0.000000,
            0.000000, 2.414213, 0.000000, 0.000000,
            0.000000, 0.000000, -1.002002, -1.000000,
            0.000000, 0.000000, -0.200200, 0.000000
        };

        template<drawable buffer_type>
        void draw(size_t size);

        void preset();
        void preset_offset();

        public:
        void * p;
        void init_shaders()
        {
            fill_shader = new Shader("shaders/fill_shader.vs",
                                 "shaders/fill_shader.fs");
            line_shader = new Shader("shaders/line_shader.vs",
                                 "shaders/line_shader.fs");
            vect_shader = new Shader("shaders/vect_shader.vs",
                                 "shaders/vect_shader.fs");
        }
        static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
        {
            float *arg;
            arg = (float *)glfwGetWindowUserPointer(window);
            if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
                glfwSetWindowShouldClose(window, true);
            if (key == GLFW_KEY_W && action == GLFW_PRESS)
                if (arg[1] < 2) arg[1] += 4;
            if (key == GLFW_KEY_S && action == GLFW_PRESS)
                if (arg[1] > -2) arg[1] -= 4;
            if (key == GLFW_KEY_A && action == GLFW_PRESS)
            {
                if (arg[0] > 0 && arg[2] == 0)
                    arg[2] += 8;
                else if (arg[0] >= 0 && arg[2] >=0)
                    arg[0] -= 6;
                else if (arg[0] < 0 && arg[2] >=0)
                    arg[2] -= 8;
                else if (arg[0] <= 0 && arg[2] <= 0)
                    arg[0] += 6;
                else if (arg[0] >= 0 && arg[2] <= 0)
                    arg[2] += 8;
            }
            if (key == GLFW_KEY_D && action == GLFW_PRESS)
            {
                if (arg[2] > 0 && arg[0] == 0)
                    arg[0] += 6;
                else if (arg[2] >= 0 && arg[0] >=0)
                    arg[2] -= 8;
                else if (arg[2] < 0 && arg[0] >=0)
                    arg[0] -= 6;
                else if (arg[2] <= 0 && arg[0] <= 0)
                    arg[2] += 8;
                else if (arg[2] >= 0 && arg[0] <= 0)
                    arg[0] += 6;
            }

            if (key == GLFW_KEY_U && action == GLFW_PRESS)
                if (arg[6] <= 1.) arg[6] += 0.1;
            if (key == GLFW_KEY_J && action == GLFW_PRESS)
                if (arg[6] >= -1.) arg[6] -= 0.1;
        }

        void draw_wrapper(size_t size, int i);
        void pre_loop();
        void post_loop();
        grid_viewer();
    };

    template<>
    void grid_viewer::draw<grid>(size_t size)
    {
        fill_shader->use();
        fill_shader->setMat4("mvp", mvp);
        fill_shader->setVec4("clip", clip);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        preset_offset();

        glDrawElements(GL_TRIANGLES, size, GL_UNSIGNED_INT, 0);

        glDisable(GL_POLYGON_OFFSET_FILL);

        line_shader->use();
        line_shader->setMat4("mvp", mvp);
        line_shader->setVec4("clip", clip);
        glLineWidth(5.0);
        glPointSize(5.0);

        preset_offset();

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        //glDrawElements(GL_LINES, size, GL_UNSIGNED_INT, 0);
        glDrawElements(GL_TRIANGLES, size, GL_UNSIGNED_INT, 0);
        return;
    }

    template<>
    void grid_viewer::draw<object>(size_t size)
    {
        /*line_shader->use();
        line_shader->setMat4("mvp", mvp);
        line_shader->setVec4("clip", clip);
        glLineWidth(5.0);
        glPointSize(5.0);

        preset();

        glDrawElements(GL_TRIANGLES, size, GL_UNSIGNED_INT, 0);
        return;*/
        fill_shader->use();
        fill_shader->setMat4("mvp", mvp);
        fill_shader->setVec4("clip", clip);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        preset_offset();

        glDrawElements(GL_TRIANGLES, size, GL_UNSIGNED_INT, 0);

        glDisable(GL_POLYGON_OFFSET_FILL);

        line_shader->use();
        line_shader->setMat4("mvp", mvp);
        line_shader->setVec4("clip", clip);
        glLineWidth(5.0);
        glPointSize(5.0);

        preset_offset();

        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

        //glDrawElements(GL_LINES, size, GL_UNSIGNED_INT, 0);
        glDrawElements(GL_TRIANGLES, size, GL_UNSIGNED_INT, 0);
        return;
    }

    template<>
    void grid_viewer::draw<vector_field>(size_t size)
    {
        vect_shader->use();
        vect_shader->setMat4("mvp", mvp);
        vect_shader->setVec4("clip", clip);

        preset();

        glDrawElements(GL_LINES, size, GL_UNSIGNED_INT, 0);
        return;
    }

    void grid_viewer::draw_wrapper(size_t size, int i)
    {
        if (i == 0)
            return draw<grid>(size);
        if (i == 1)
            return draw<object>(size);
        if (i == 2)
            return draw<vector_field>(size);
    }

    void grid_viewer::pre_loop()
    {
        glClearColor(0.7f, 0.7f, 0.7f, 0.7f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 projection;
        memcpy( glm::value_ptr(projection), aaa, sizeof(aaa));

        glm::mat4 view = glm::lookAt(
            camera_pos,
            glm::vec3(0,0,0),
            glm::vec3(0,1,0)
        );

        mvp = projection * view;
        //glEnable(GL_CLIP_PLANE0);
        return;
    }

    void grid_viewer::post_loop()
    {
        //glDisable(GL_CLIP_PLANE0);
        return;
    }

    void grid_viewer::preset()
    {
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
                              3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        return;
    }

    void grid_viewer::preset_offset()
    {
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
                              6 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE,
                              6 * sizeof(float),
                              (void*)(3*sizeof(float)) );
        glEnableVertexAttribArray(1);
        return;
    }

    grid_viewer::grid_viewer()
    {
        camera_pos = glm::vec3{-6., 4., 8.};
        clip = glm::vec4{0., 0., -1., 2.};
        p = (void*) &camera_pos;
    }
}
