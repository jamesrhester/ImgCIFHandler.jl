# These macros will become more sophisticated
# with time to make a nice printout. For now
# they simply collect the tests into three
# lists.
macro noimgcheck(description, check)
    # extract name
    #func_name = expr.args[1].args[1] #:=->:call
    #push!(test_list,func_name)
    quote
        x = $(esc(check))
        push!(test_list,($(esc(description)),x))
    end
end

macro imgcheck(description, check)
    quote
        x = $(esc(check))
        push!(test_list_with_img,($(esc(description)),x))
    end
end
        
macro fullcheck(description, check)
    quote
        x = $(esc(check))
        push!(test_full_list,($(esc(description)), x))
    end
end

#macro plaincheck(testname,testexpr,message)
#end
