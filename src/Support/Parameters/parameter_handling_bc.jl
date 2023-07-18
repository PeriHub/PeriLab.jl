function get_bcs(params)

    check = check_element(params["Boundary Conditions"])
    if check
        for entry in params["Boundary Conditions"]
            if ".txt" in params["Boundary Conditions"][entry]
            end
        end
    end
end