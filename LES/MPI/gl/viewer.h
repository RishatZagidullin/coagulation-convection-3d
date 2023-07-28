#pragma once
#include <GLFW/glfw3.h>

namespace gl
{
    template<int n_objects, typename rendering>
    class viewer
    {
        private:
        short SCR_WIDTH;
        short SCR_HEIGHT;
        unsigned int VAO;
        rendering r;
        public:
        GLFWwindow *window;
        unsigned int buffer_objects [n_objects*2];

        template<int buffer_type, typename T> 
        void buffer(size_t const &size, T const * obj = nullptr, size_t const &i = 0)
        {
            glBindBuffer(buffer_type, buffer_objects[i]);
            glBufferData(buffer_type, sizeof(T)*size, obj, GL_DYNAMIC_DRAW);
            return;
        }

        void view(int const sizes [n_objects], int time = -1);

        viewer(short SCR_WIDTH, short SCR_HEIGHT, bool hide_window = false);
        ~viewer();
    };

    template<int n_objects, typename rendering>
    void viewer<n_objects, rendering>::view(int const sizes [n_objects], int time)
    {
        r.pre_loop();
        for (int i = 0; i < n_objects; i++)
        {
            glBindBuffer(GL_ARRAY_BUFFER, buffer_objects[2*i]);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_objects[2*i+1]);

            r.draw_wrapper(sizes[i], i);
        }
        r.post_loop();
        glfwSwapBuffers(window);
        glfwPollEvents();

        if (time != -1)
        {
            int * buff = new int [SCR_WIDTH * SCR_HEIGHT * 3];
            glReadPixels(0, 0, SCR_WIDTH, SCR_HEIGHT, GL_BGR,
                         GL_UNSIGNED_BYTE, buff);
            std::string filename = std::string("imgs/data_")
                    + std::to_string(time) + std::string(".tga");
            FILE *out = std::fopen(filename.c_str(), "w");
            short head[] {0,2,0,0,0,0,SCR_WIDTH,SCR_HEIGHT,24};
            fwrite(&head, sizeof(head), 1, out);
            fwrite(buff, 3*SCR_WIDTH*SCR_HEIGHT, 1, out);
            fclose(out);
            delete [] buff;
        }
        return;
    }

    template<int n_objects, typename rendering>
    viewer<n_objects, rendering>::viewer(short SCR_WIDTH, short SCR_HEIGHT, bool hide_window)
    {
        this->SCR_WIDTH = SCR_WIDTH;
        this->SCR_HEIGHT = SCR_HEIGHT;
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        if (hide_window)
            glfwWindowHint(GLFW_VISIBLE, GL_FALSE);

        window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT,
                                  "", NULL, NULL);
        glfwMakeContextCurrent(window);
        gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
        const GLubyte* vendor = glGetString(GL_VENDOR);
        const GLubyte* renderer = glGetString(GL_RENDERER);
        std::cout << vendor << "  " << renderer << "\n";

        r.init_shaders();

        glEnable(GL_DEPTH_TEST);
        glGenVertexArrays(1, &VAO);
        
        for (int i = 0; i < n_objects; i++)
        {
            glGenBuffers(1, &buffer_objects[2*i]);
            glGenBuffers(1, &buffer_objects[2*i+1]);
        }
        glBindVertexArray(VAO);

        glfwSetWindowUserPointer(window, (void *)r.p);
        glfwSetKeyCallback(window, r.key_callback);
        //if (glGetError() == GL_NO_ERROR)
        //    std::cout << "good\n";
    }

    template<int n_objects, typename rendering>
    viewer<n_objects, rendering>::~viewer()
    {
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDeleteVertexArrays(1, &VAO);
        for (int i = 0; i < n_objects; i++)
        {
            glDeleteBuffers(1, &buffer_objects[2*i]);
            glDeleteBuffers(1, &buffer_objects[2*i+1]);
        }
        glfwTerminate();
    }
}
